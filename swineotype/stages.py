from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from swineotype.blast import run_blast, make_db_if_needed

from swineotype.utils import ensure_tool, gzip_file

# Stage-2 needs the aligned sequences so the diagnostic base can be read out of
# the alignment itself rather than by arithmetic on start coordinates.
RESOLVER_OUTFMT = ("6 qseqid sseqid pident length qlen evalue bitscore "
                   "qstart qend sstart send qseq sseq")


def reverse_complement(base: str) -> str:
    """Returns the reverse complement of a single DNA base."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return comp.get(base.upper(), 'N')

def parse_whitelist_headers(fasta_path: str):
    allele_to_type = {}
    allele_to_geneclass = {}
    with open(fasta_path, "r") as fh:
        for line in fh:
            if not line.startswith(">"): continue
            h = line[1:].strip()
            allele_id = h.split()[0]
            st = None
            for tok in h.split():
                if tok.startswith("[type_id=") and tok.endswith("]"):
                    st = tok[len("[type_id="):-1]; break
            allele_to_type[allele_id] = st
            low = allele_id.lower()
            geneclass = "wzy" if "wzy" in low else ("wzx" if "wzx" in low else None)
            allele_to_geneclass[allele_id] = geneclass
    return allele_to_type, allele_to_geneclass

def stage1_score(assembly_fa: str, whitelist_fa: str, threads: int, run_dir: Path, config: dict):
    ensure_tool("blastn"); ensure_tool("makeblastdb")
    allele_to_type, allele_to_geneclass = parse_whitelist_headers(whitelist_fa)
    db_prefix = make_db_if_needed(assembly_fa, config["tmp_dir"])
    outfmt = "6 qseqid sseqid pident length qlen evalue bitscore qstart qend sstart send"
    tsv_text = run_blast(whitelist_fa, db_prefix, threads, outfmt)
    stage1_tsv = run_dir / "wzxwzy_vs_asm.tsv"
    if config["keep_debug"]:
        stage1_tsv.write_text(tsv_text + ("\n" if tsv_text else ""))

        if config["gzip_debug"]:
            gzip_file(stage1_tsv)

    hits = defaultdict(list)
    for line in filter(None, tsv_text.splitlines()):
        parts = line.split("\t")
        if len(parts) < 11: continue
        qseqid, sseqid, pident, length, qlen, evalue, bitscore, qstart, qend, sstart, send = parts
        # Group by qseqid ONLY (ignore sseqid/contig to handle genes split across contigs)
        hits[qseqid].append({
            "pident": float(pident),
            "length": int(length),
            "qlen": int(qlen),
            "bitscore": float(bitscore),
            "qstart": int(qstart),
            "qend": int(qend)
        })

    score_by_type = defaultdict(float)

    best_rejected_info = None

    for qseqid, hsp_list in hits.items():
        # 1. Calculate total query coverage by merging intervals
        intervals = sorted([(h["qstart"], h["qend"]) for h in hsp_list])
        merged = []
        if intervals:
            curr_start, curr_end = intervals[0]
            for next_start, next_end in intervals[1:]:
                if next_start <= curr_end + 1: # overlapping or adjacent
                    curr_end = max(curr_end, next_end)
                else:
                    merged.append((curr_start, curr_end))
                    curr_start, curr_end = next_start, next_end
            merged.append((curr_start, curr_end))

        covered_len = sum(e - s + 1 for s, e in merged)
        qlen = hsp_list[0]["qlen"]
        coverage = (covered_len / qlen) if qlen else 0.0

        # 2. Calculate weighted PID (length-weighted)
        total_aligned_len = sum(h["length"] for h in hsp_list)
        weighted_pid_sum = sum(h["pident"] * h["length"] for h in hsp_list)
        avg_pid = (weighted_pid_sum / total_aligned_len) if total_aligned_len else 0.0

        # 3. Filter
        if avg_pid < config["min_pid"] or coverage < config["min_cov"]:
            # Diagnostic for best rejected
            info = f"{qseqid}: Cov={coverage:.2f} Pid={avg_pid:.1f}"
            if best_rejected_info is None or coverage > float(best_rejected_info.split("Cov=")[1].split()[0]):
                 best_rejected_info = info
            continue

        # 4. Score this allele: sum bitscores over a greedy, non-overlapping
        # set of HSPs (best-scoring first).
        #
        # Summing over ALL HSPs let a duplicated or repeated copy of a gene in
        # the assembly count twice, inflating that serotype's score enough to
        # push it to the top of the ranking. Taking only the single best HSP
        # would fix that but would under-score a gene legitimately split across
        # two contigs. Requiring the counted HSPs not to overlap in QUERY space
        # handles both: split genes contribute each of their parts, duplicate
        # copies of the same region do not.
        allele_bitscore = 0.0
        counted: list[tuple[int, int]] = []
        for h in sorted(hsp_list, key=lambda x: -x["bitscore"]):
            s, e = min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])
            if any(s <= c_e and e >= c_s for c_s, c_e in counted):
                continue
            counted.append((s, e))
            allele_bitscore += h["bitscore"]

        st = allele_to_type.get(qseqid)
        if st:
            score_by_type[st] += allele_bitscore

    if not score_by_type and best_rejected_info:
        print(f"[DEBUG] No hits passed filter. Best rejected: {best_rejected_info}")


    ordered = sorted(score_by_type.items(), key=lambda kv: kv[1], reverse=True)
    top, top_score = (ordered[0][0], ordered[0][1]) if ordered else (None, 0.0)
    second, second_score = (ordered[1][0], ordered[1][1]) if len(ordered) > 1 else (None, 0.0)
    total = sum(score_by_type.values())
    fraction, delta = (top_score/total if total else 0.0), top_score-second_score
    decisive = (fraction >= config["plurality"]) and (delta >= config["delta"])
    must_stage2_for_pair = bool(top in config["ambig_set"])
    return {"scores":score_by_type,"top":top,"second":second,"fraction":fraction,
            "delta":delta,"decisive":decisive,"must_stage2_for_pair":must_stage2_for_pair}


def base_at_query_pos(qseq: str, sseq: str, qstart: int, sstart: int, send: int, pos: int):
    """Read the subject base aligned to query position ``pos``.

    Walks the gapped alignment column by column instead of computing
    ``sstart + (pos - qstart)``. That arithmetic silently assumed an ungapped
    HSP, so a single indel anywhere between the HSP start and the diagnostic
    site shifted the read-out by one base -- the characteristic ONT error mode,
    and enough to flip a serotype call.

    BLAST reports qseq/sseq in alignment orientation, so for a minus-strand
    subject hit sseq is already reverse-complemented onto the query strand and
    needs no further complementing.

    Returns ``(base, contig_pos, strand)``. ``base`` is ``"-"`` when the
    subject carries a deletion at the diagnostic site, or ``None`` if the
    alignment does not actually reach ``pos``.
    """
    strand = "+" if sstart <= send else "-"
    sstep = 1 if strand == "+" else -1

    qpos, spos = qstart, sstart
    for qchar, schar in zip(qseq, sseq):
        q_gap = (qchar == "-")
        s_gap = (schar == "-")
        if not q_gap and qpos == pos:
            return ("-" if s_gap else schar.upper()), (None if s_gap else spos), strand
        if not q_gap:
            qpos += 1
        if not s_gap:
            spos += sstep
    return None, None, strand


def stage2_resolver_call(assembly_fa: str, resolver_refs_fa: str, threads: int, run_dir: Path, config: dict, allowed_pair: str|None=None):
    ensure_tool("blastn"); ensure_tool("makeblastdb")
    db_prefix = make_db_if_needed(assembly_fa, config["tmp_dir"])
    tsv_text = run_blast(resolver_refs_fa, db_prefix, threads, RESOLVER_OUTFMT)

    if config["keep_debug"]:
        stage2_tsv = run_dir / "resolver_vs_asm.tsv"
        stage2_tsv.write_text(tsv_text + ("\n" if tsv_text else ""))
        if config["gzip_debug"]:
            gzip_file(stage2_tsv)


    best = None
    for line in filter(None, tsv_text.splitlines()):
        parts = line.split("\t")
        if len(parts) < 13: continue
        (qseqid, sseqid, pident, length, qlen, evalue, bitscore,
         qstart, qend, sstart, send, qseq, sseq) = parts[:13]
        pident, length, qstart, qend, sstart, send, bitscore = float(pident), int(length), int(qstart), int(qend), int(sstart), int(send), float(bitscore)
        meta = parse_resolver_meta(qseqid)
        if allowed_pair and meta["pair"] != allowed_pair: continue
        pos = meta["pos"]
        spans = (qstart <= pos <= qend) or (qend <= pos <= qstart)
        if not spans: continue
        if pident < config["min_res_pid"] or length < config["min_res_alen"]: continue
        base, tpos, strand = base_at_query_pos(qseq, sseq, qstart, sstart, send, pos)
        if base is None: continue
        ev = {"ref_id":qseqid,"contig":sseqid,"contig_pos":tpos,"strand":strand,
              "pident":pident,"length":length,"bitscore":bitscore,"pair":meta["pair"],"base":base}
        if best is None or bitscore > best[0]: best = (bitscore, ev)
    if not best: return None
    return best[1]


def interpret_resolver(ev: dict|None, config: dict) -> str|None:
    if ev is None: return None
    meta = parse_resolver_meta(ev["ref_id"])
    base = ev["base"]
    if base == "G":
        return meta.get("G_serotype")
    elif base in ("C", "T"):
        return meta.get("CT_serotype")
    return None

def parse_resolver_meta(qid: str):
    """
    Parse resolver FASTA IDs of the form:
      >id|pair=1_vs_14|pos=492|G_serotype=14|CT_serotype=1
    Returns dict with keys: pair, pos, G_serotype, CT_serotype

    `pos` is 1-based and relative to THAT reference, which is why the 1/14 and
    2/(1/2) references declare different positions (492 vs 483) for the same
    CpsK residue -- the 1/14 references carry 9 extra bases at the 5' end.
    """
    meta = {"pair": None, "pos": None, "G_serotype": None, "CT_serotype": None}
    for tok in qid.split("|")[1:]:
        if tok.startswith("pair="):
            meta["pair"] = tok.split("=", 1)[1]
        elif tok.startswith("pos="):
            meta["pos"] = int(tok.split("=", 1)[1])
        elif tok.startswith("G_serotype="):
            meta["G_serotype"] = tok.split("=", 1)[1]
        elif tok.startswith("CT_serotype="):
            meta["CT_serotype"] = tok.split("=", 1)[1]
    return meta
