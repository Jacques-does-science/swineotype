from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from swineotype.blast import run_blast, make_db_if_needed, run

from swineotype.utils import ensure_tool, gzip_file


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

    score_by_type = defaultdict(float)
    for line in filter(None, tsv_text.splitlines()):
        qseqid, sseqid, pident, length, qlen, *_rest = line.split("\t")
        pident, length, qlen = float(pident), int(length), int(qlen)
        coverage = (length / qlen) if qlen else 0
        if pident < config["min_pid"] or coverage < config["min_cov"]: continue
        bitscore = float(_rest[1])
        st = allele_to_type.get(qseqid)
        if st: score_by_type[st] += bitscore

    ordered = sorted(score_by_type.items(), key=lambda kv: kv[1], reverse=True)
    top, top_score = (ordered[0][0], ordered[0][1]) if ordered else (None, 0.0)
    second, second_score = (ordered[1][0], ordered[1][1]) if len(ordered) > 1 else (None, 0.0)
    total = sum(score_by_type.values())
    fraction, delta = (top_score/total if total else 0.0), top_score-second_score
    decisive = (fraction >= config["plurality"]) and (delta >= config["delta"])
    must_stage2_for_pair = bool(top in config["ambig_set"])
    return {"scores":score_by_type,"top":top,"second":second,"fraction":fraction,
            "delta":delta,"decisive":decisive,"must_stage2_for_pair":must_stage2_for_pair}


def stage2_resolver_call(assembly_fa: str, resolver_refs_fa: str, threads: int, run_dir: Path, config: dict, allowed_pair: str|None=None):
    ensure_tool("blastn"); ensure_tool("makeblastdb"); ensure_tool("samtools")
    db_prefix = make_db_if_needed(assembly_fa, config["tmp_dir"])
    outfmt = "6 qseqid sseqid pident length qlen evalue bitscore qstart qend sstart send"
    tsv_text = run_blast(resolver_refs_fa, db_prefix, threads, outfmt)

    if config["keep_debug"]:
        stage2_tsv = run_dir / "resolver_vs_asm.tsv"
        stage2_tsv.write_text(tsv_text + ("\n" if tsv_text else ""))
        if config["gzip_debug"]:
            gzip_file(stage2_tsv)


    best = None
    for line in filter(None, tsv_text.splitlines()):
        qseqid, sseqid, pident, length, qlen, evalue, bitscore, qstart, qend, sstart, send = line.split("\t")
        pident, length, qstart, qend, sstart, send, bitscore = float(pident), int(length), int(qstart), int(qend), int(sstart), int(send), float(bitscore)
        meta = parse_resolver_meta(qseqid)
        if allowed_pair and meta["pair"] != allowed_pair: continue
        pos = meta["pos"]
        spans = (qstart <= pos <= qend) or (qend <= pos <= qstart)
        if not spans: continue
        if pident < config["min_res_pid"] or length < config["min_res_alen"]: continue
        qoff = pos - qstart
        if sstart <= send:
            strand = "+"; tpos = sstart + qoff
        else:
            strand = "-"; tpos = sstart - qoff
        ev = {"ref_id":qseqid,"contig":sseqid,"contig_pos":tpos,"strand":strand,
              "pident":pident,"length":length,"bitscore":bitscore,"pair":meta["pair"],"base":None}
        if best is None or bitscore > best[0]: best = (bitscore, ev)
    if not best: return None
    ev = best[1]
    region = f"{ev['contig']}:{ev['contig_pos']}-{ev['contig_pos']}"
    fa = run(["samtools","faidx",assembly_fa,region])
    lines = [ln.strip() for ln in fa.splitlines()]
    ev["base"] = lines[1].strip().upper() if len(lines)>1 else "N"
    return ev


def interpret_resolver(ev: dict|None, config: dict) -> str|None:
    if ev is None: return None
    meta = parse_resolver_meta(ev["ref_id"])
    if meta["pair"] == "1_vs_14":
        return meta["B"] if ev["base"] == meta["baseA"] else meta["A"]
    if meta["pair"] == "2_vs_1_2":
        return meta["B"] if ev["base"] == meta["baseA"] else meta["A"]
    return None

def parse_resolver_meta(qid: str):
    """
    Parse resolver FASTA IDs of the form:
      >id|pair=1_vs_14|pos=483|baseA=G|A=1|B=14
    Returns dict with keys: pair, pos, baseA, A, B
    """
    meta = {"pair": None, "pos": None, "baseA": "G", "A": None, "B": None}
    for tok in qid.split("|")[1:]:
        if tok.startswith("pair="):   meta["pair"]   = tok.split("=",1)[1]
        elif tok.startswith("pos="):  meta["pos"]    = int(tok.split("=",1)[1])
        elif tok.startswith("baseA="):meta["baseA"]  = tok.split("=",1)[1].upper()
        elif tok.startswith("A="):    meta["A"]      = tok.split("=",1)[1]
        elif tok.startswith("B="):    meta["B"]      = tok.split("=",1)[1]
    return meta
