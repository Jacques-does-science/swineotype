from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from swineotype.blast import run_blast, make_db_if_needed

from swineotype.utils import ensure_tool, gzip_file

# Stage-2 needs the aligned sequences so the diagnostic codon can be read out
# of the alignment itself rather than by arithmetic on start coordinates.
RESOLVER_OUTFMT = ("6 qseqid sseqid pident length qlen evalue bitscore "
                   "qstart qend sstart send qseq sseq")

STAGE1_OUTFMT = "6 qseqid sseqid pident length qlen evalue bitscore qstart qend sstart send"

DNA = frozenset("ACGT")
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

# The only diagnostic codons the scheme defines. Trp161 (TGG) -> Gal -> the
# "G" serotype of the pair; Cys161 (TGT/TGC) -> GalNAc -> the "CT" serotype.
# Anything else -- AGG, CGG, TAG, an ambiguity code -- is a state the scheme
# has no interpretation for and must not be coerced into one.
TRP_CODONS = frozenset({"TGG"})
CYS_CODONS = frozenset({"TGT", "TGC"})
RESOLVABLE_CODONS = TRP_CODONS | CYS_CODONS

# Triplet read-out states.
TRIPLET_OK = "OK"                        # complete, contiguous, TGG/TGC/TGT
TRIPLET_UNEXPECTED = "UNEXPECTED_CODON"  # complete and contiguous, but not one of those
TRIPLET_AMBIGUOUS = "AMBIGUOUS"          # contains a non-ACGT character
TRIPLET_DELETED = "DELETED"              # subject carries a deletion in the codon
TRIPLET_INCOMPLETE = "INCOMPLETE"        # alignment does not reach all three positions
TRIPLET_NONCONTIGUOUS = "NON_CONTIGUOUS" # subject insertion inside the codon

# Coding-sequence integrity states. UNASSESSED is not a synonym for intact:
# it means the alignment did not carry enough of the gene to judge.
CODING_INTACT = "INTACT"
CODING_DISRUPTED = "DISRUPTED"
CODING_UNASSESSED = "UNASSESSED"

# Stage-2 outcomes.
S2_SKIPPED = "SKIPPED"
S2_NO_HSP = "NO_HSP_OR_LOW_QUAL"
S2_OK = "OK"
S2_CONFLICT = "CONFLICTING_COPIES"
S2_INVALID = "INVALID_TRIPLET"
S2_CODING_DISRUPTED = "CODING_DISRUPTED"


# There is deliberately no reverse_complement() here. BLAST reports sseq in
# alignment orientation, already on the query strand, so a minus-strand hit
# needs no complementing -- only its subject COORDINATE runs backwards. A
# helper of that name sitting next to this code invites exactly the
# double-complement bug the read-out is built to avoid.


def header_tag(header: str, key: str) -> str | None:
    """Value of a `[key=value]` tag in a FASTA header, or None.

    Values may contain spaces (species binomials, literature references), so
    this cannot be a whitespace split.
    """
    match = re.search(r"\[" + re.escape(key) + r"=([^\]]*)\]", header)
    return match.group(1) if match else None


def parse_whitelist_headers(fasta_path: str):
    """Returns (allele -> cps type, allele -> wzx/wzy, cps type -> species).

    The species comes from each reference's own `[species=...]` tag rather than
    being assumed. Former S. suis serotypes 20, 22 and 26 are Streptococcus
    parasuis, 33 is S. ruminantium, and 32 and 34 are S. orisratti; their cps
    loci are still in the panel, so a hit to one identifies a different
    organism, not an S. suis serotype.
    """
    allele_to_type = {}
    allele_to_geneclass = {}
    type_to_species = {}
    with open(fasta_path, "r") as fh:
        for line in fh:
            if not line.startswith(">"): continue
            h = line[1:].strip()
            allele_id = h.split()[0]
            st = header_tag(h, "type_id")
            allele_to_type[allele_id] = st
            if st is not None:
                type_to_species[st] = header_tag(h, "species")
            low = allele_id.lower()
            geneclass = "wzy" if "wzy" in low else ("wzx" if "wzx" in low else None)
            allele_to_geneclass[allele_id] = geneclass
    return allele_to_type, allele_to_geneclass, type_to_species


# ---------------------------------------------------------------------------
# Stage 1: HSP aggregation
#
# Coverage, identity and score are all computed from ONE selected set of HSPs
# -- a single physical copy of the gene, or an explicitly-flagged chain across
# contigs. The previous code computed each of the three from a different set:
# coverage was unioned over every HSP, identity was length-weighted over every
# HSP including overlapping copies, and only the score used a greedy
# non-overlapping subset. A full-length 100%-identity match therefore passed on
# its own but failed once five overlapping 80%-identity copies were added,
# because they dragged the reported identity below the threshold while
# contributing nothing to the score.
# ---------------------------------------------------------------------------

def _q_span(h) -> tuple[int, int]:
    return min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])


def _q_overlap(a, b) -> bool:
    a_lo, a_hi = _q_span(a)
    b_lo, b_hi = _q_span(b)
    return a_lo <= b_hi and b_lo <= a_hi


def _collinear(a, b) -> bool:
    """Do two HSPs sit in the same order in query and in subject space?

    Two fragments of one gene advance together along query and subject (in the
    strand's direction). Two fragments in inconsistent order are a rearranged
    or repeated region, not one copy.
    """
    dq = a["qstart"] - b["qstart"]
    ds = a["sstart"] - b["sstart"]
    if dq == 0 or ds == 0:
        return False
    same_direction = (dq > 0) == (ds > 0)
    return same_direction if a["strand"] == "+" else not same_direction


def _merge_intervals(intervals):
    merged = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def _metrics(hsps, qlen) -> dict:
    """Coverage, identity and score for ONE selected set of HSPs."""
    covered = sum(e - s + 1 for s, e in _merge_intervals([_q_span(h) for h in hsps]))
    aligned = sum(h["length"] for h in hsps)
    return {
        "coverage": (covered / qlen) if qlen else 0.0,
        "identity": (sum(h["pident"] * h["length"] for h in hsps) / aligned) if aligned else 0.0,
        "score": sum(h["bitscore"] for h in hsps),
        "aligned_len": aligned,
        "covered_len": covered,
        "qlen": qlen,
    }


def group_into_copies(hsps: list[dict]) -> list[list[dict]]:
    """Partition one query's HSPs into distinct physical copies.

    Two HSPs belong to the same copy when they are on the same contig, on the
    same strand, do not overlap in query space, and are collinear with every
    HSP already in that copy. Best-scoring HSP first, so the strongest copy is
    seeded first and a weaker overlapping alternative is pushed into a copy of
    its own instead of being blended into the strong one.
    """
    copies: list[list[dict]] = []
    for h in sorted(hsps, key=lambda x: (-x["bitscore"], _q_span(x))):
        for members in copies:
            m0 = members[0]
            if m0["sseqid"] != h["sseqid"] or m0["strand"] != h["strand"]:
                continue
            if any(_q_overlap(h, m) for m in members):
                continue
            if not all(_collinear(h, m) for m in members):
                continue
            members.append(h)
            break
        else:
            copies.append([h])
    return copies


def allele_evidence(hsps: list[dict], qlen: int, min_cov: float) -> dict:
    """The evidence one reference allele contributes, as a single copy.

    Policy ("best consistent copy, then an explicit split chain"):

    1. Partition the HSPs into physical copies (see group_into_copies).
    2. The allele's evidence is the highest-scoring copy. Coverage, identity
       and score all come from that copy and no other, so a weaker overlapping
       duplicate can neither inflate the score nor dilute the identity.
    3. Only if that copy alone does not reach ``min_cov`` is a chain across
       copies attempted, adding HSPs that cover disjoint query intervals,
       best-scoring first. That keeps genuinely split-contig genes callable,
       but the result is flagged ``split=True`` and carries every contig it
       drew on, because fragments from unrelated locations do not by
       themselves prove one intact gene.
    """
    copies = group_into_copies(hsps)
    scored = sorted(((_metrics(c, qlen), c) for c in copies),
                    key=lambda mc: -mc[0]["score"])
    best_metrics, best_hsps = scored[0]
    chosen, split = list(best_hsps), False

    if best_metrics["coverage"] < min_cov and len(copies) > 1:
        others = [h for _, c in scored[1:] for h in c]
        for h in sorted(others, key=lambda x: -x["bitscore"]):
            if not any(_q_overlap(h, c) for c in chosen):
                chosen.append(h)
        candidate = _metrics(chosen, qlen)
        if candidate["coverage"] > best_metrics["coverage"]:
            best_metrics, split = candidate, True
        else:
            chosen = list(best_hsps)

    contigs = sorted({h["sseqid"] for h in chosen})
    return {
        **best_metrics,
        "split": split,
        "contigs": contigs,
        "n_copies": len(copies),
        "strands": sorted({h["strand"] for h in chosen}),
        "n_hsps": len(chosen),
        "subject_spans": sorted((h["sseqid"], min(h["sstart"], h["send"]),
                                 max(h["sstart"], h["send"])) for h in chosen),
    }


def family_of(serotype: str | None, config: dict) -> str | None:
    """The evidence family a cps type belongs to.

    {1,14} and {2,1/2} each share a locus and differ only at the cpsK
    diagnostic site, so their references are near-duplicates of one another and
    must not be counted as independent support. Every other type is its own
    family.
    """
    if serotype is None:
        return None
    if serotype in config.get("pair_1_14", ()):
        return "1_vs_14"
    if serotype in config.get("pair_2_1_2", ()):
        return "2_vs_1_2"
    return f"type:{serotype}"


RESOLVABLE_FAMILIES = ("1_vs_14", "2_vs_1_2")

FAMILY_MEMBERS = {"1_vs_14": ("1", "14"), "2_vs_1_2": ("2", "1/2")}


def family_members(family: str | None, config: dict | None = None) -> tuple[str, ...]:
    """The cps types a family contains, honouring a reconfigured pair set."""
    if family is None:
        return ()
    if family.startswith("type:"):
        return (family.split(":", 1)[1],)
    key = {"1_vs_14": "pair_1_14", "2_vs_1_2": "pair_2_1_2"}.get(family)
    configured = (config or {}).get(key) if key else None
    if configured:
        return tuple(sorted(configured, key=lambda s: (len(s), s)))
    return FAMILY_MEMBERS.get(family, ())


def family_label(family: str | None, config: dict | None = None) -> str:
    """Human-readable family-level result, e.g. `1 or 14`."""
    members = family_members(family, config)
    if members:
        return " or ".join(members)
    return family or ""


def stage1_score(assembly_fa: str, whitelist_fa: str, threads: int, run_dir: Path, config: dict):
    ensure_tool("blastn"); ensure_tool("makeblastdb")
    allele_to_type, allele_to_geneclass, type_to_species = parse_whitelist_headers(whitelist_fa)
    db_prefix = make_db_if_needed(assembly_fa, config["tmp_dir"])
    tsv_text = run_blast(whitelist_fa, db_prefix, threads, STAGE1_OUTFMT)
    stage1_tsv = run_dir / "wzxwzy_vs_asm.tsv"
    if config["keep_debug"]:
        stage1_tsv.write_text(tsv_text + ("\n" if tsv_text else ""))

        if config["gzip_debug"]:
            gzip_file(stage1_tsv)

    hits = defaultdict(list)
    for line in filter(None, tsv_text.splitlines()):
        parts = line.split("\t")
        if len(parts) < 11: continue
        qseqid, sseqid, pident, length, qlen, evalue, bitscore, qstart, qend, sstart, send = parts[:11]
        sstart, send = int(sstart), int(send)
        # Group by qseqid, but keep the subject coordinates and orientation so
        # copy identity survives into the aggregation. Grouping by qseqid ONLY
        # and discarding sseqid -- which is what the tool used to do -- made a
        # duplicated gene indistinguishable from a gene split across contigs.
        hits[qseqid].append({
            "pident": float(pident),
            "length": int(length),
            "qlen": int(qlen),
            "bitscore": float(bitscore),
            "qstart": int(qstart),
            "qend": int(qend),
            "sseqid": sseqid,
            "sstart": sstart,
            "send": send,
            "strand": "+" if sstart <= send else "-",
        })

    evidence_by_allele: dict[str, dict] = {}
    best_rejected = None

    for qseqid, hsp_list in hits.items():
        ev = allele_evidence(hsp_list, hsp_list[0]["qlen"], config["min_cov"])
        if ev["identity"] < config["min_pid"] or ev["coverage"] < config["min_cov"]:
            info = f"{qseqid}: Cov={ev['coverage']:.2f} Pid={ev['identity']:.1f}"
            if best_rejected is None or ev["coverage"] > best_rejected[0]:
                best_rejected = (ev["coverage"], info)
            continue
        evidence_by_allele[qseqid] = ev

    # Best score per (type, gene class) and per (family, gene class). Taking the
    # best rather than the sum is what stops near-duplicate references inside
    # one family -- serotype 1's and serotype 14's wzy, which differ at a single
    # site -- from being counted as two independent pieces of evidence.
    type_best: dict[str, dict[str, float]] = defaultdict(dict)
    family_best: dict[str, dict[str, float]] = defaultdict(dict)
    genes_by_type = defaultdict(set)
    family_types = defaultdict(set)
    family_evidence = defaultdict(dict)

    for qseqid, ev in evidence_by_allele.items():
        st = allele_to_type.get(qseqid)
        if not st:
            continue
        gc = allele_to_geneclass.get(qseqid) or "other"
        fam = family_of(st, config)
        type_best[st][gc] = max(type_best[st].get(gc, 0.0), ev["score"])
        if ev["score"] >= family_best[fam].get(gc, 0.0):
            family_best[fam][gc] = ev["score"]
            family_evidence[fam][gc] = {"allele": qseqid, "type": st, **ev}
        family_types[fam].add(st)
        if allele_to_geneclass.get(qseqid):
            genes_by_type[st].add(gc)

    score_by_type = {t: sum(g.values()) for t, g in type_best.items()}
    score_by_family = {f: sum(g.values()) for f, g in family_best.items()}

    if not score_by_type and best_rejected:
        print(f"[DEBUG] No hits passed filter. Best rejected: {best_rejected[1]}")

    ranked_all = sorted(score_by_type.items(), key=lambda kv: kv[1], reverse=True)

    # Only a type whose *wzy* survived the filters may be called.
    #
    # wzx (flippase) is conserved across serotypes -- a single wzx routinely
    # matches four different references at 95-98% -- whereas wzy (polymerase)
    # is the serotype-specific gene, which is why the published multiplex PCR
    # schemes target it. Requiring wzx as well would be wrong: serotype 14 has
    # no wzx reference.
    require_wzy = bool(config.get("require_wzy"))
    eligible = [(t, s) for t, s in ranked_all if "wzy" in genes_by_type[t]] \
        if require_wzy else ranked_all

    top, top_score = eligible[0] if eligible else (None, 0.0)
    second, second_score = eligible[1] if len(eligible) > 1 else (None, 0.0)
    total = sum(s for _, s in eligible)
    fraction, delta = (top_score / total if total else 0.0), top_score - second_score
    decisive = (fraction >= config["plurality"]) and (delta >= config["delta"])

    # --- family-level confidence --------------------------------------
    #
    # Confidence between families is established BEFORE any within-family
    # resolution is attempted. Ambiguity inside {1,14} or inside {2,1/2} is
    # expected and is exactly what Stage 2 exists to settle; ambiguity between
    # a resolvable family and an unrelated serotype is not, and must not be
    # settled by reading a cpsK site that only distinguishes members of the
    # family.
    eligible_families = [
        (f, s) for f, s in sorted(score_by_family.items(), key=lambda kv: -kv[1])
        if (not require_wzy) or "wzy" in family_best[f]
    ]
    fam_top, fam_top_score = eligible_families[0] if eligible_families else (None, 0.0)
    fam_second, fam_second_score = eligible_families[1] if len(eligible_families) > 1 else (None, 0.0)
    fam_total = sum(s for _, s in eligible_families)
    fam_fraction = (fam_top_score / fam_total) if fam_total else 0.0
    fam_delta = fam_top_score - fam_second_score
    # Same two heuristics as the type-level gate, applied to families. They are
    # heuristics, not calibrated probabilities: `plurality` is the share of
    # total family evidence the leader must hold, and `delta` is the raw
    # bit-score margin it must hold over the runner-up family. Both must pass.
    fam_decisive = (fam_fraction >= config["plurality"]) and (fam_delta >= config["delta"])

    must_stage2_for_pair = fam_top in RESOLVABLE_FAMILIES

    # Best hit that is NOT callable for want of a wzy: the diagnostic for a
    # capsular locus outside the panel.
    wzx_only = [(t, s) for t, s in ranked_all if "wzy" not in genes_by_type[t]]

    return {"scores": score_by_type, "top": top, "second": second, "fraction": fraction,
            "delta": delta, "decisive": decisive,
            "must_stage2_for_pair": must_stage2_for_pair,
            "type_to_species": type_to_species,
            "top_species": type_to_species.get(top) if top else None,
            "genes_by_type": {t: sorted(g) for t, g in genes_by_type.items()},
            "top_wzx_only": wzx_only[0][0] if wzx_only else None,
            "family_scores": score_by_family,
            "family_ranked": eligible_families,
            "family_top": fam_top, "family_second": fam_second,
            "family_fraction": fam_fraction, "family_delta": fam_delta,
            "family_decisive": fam_decisive,
            "family_types": {f: sorted(t) for f, t in family_types.items()},
            "family_evidence": {f: dict(g) for f, g in family_evidence.items()},
            "allele_evidence": evidence_by_allele}


# ---------------------------------------------------------------------------
# Stage 2: gap-aware read-out of the diagnostic codon
# ---------------------------------------------------------------------------

def alignment_columns(qseq: str, sseq: str, qstart: int, sstart: int, send: int):
    """Walk a gapped HSP, yielding (qpos, spos, qchar, schar) per column.

    ``qpos``/``spos`` are None in the column where that sequence carries a gap.
    BLAST reports qseq/sseq in alignment orientation, so for a minus-strand
    subject hit sseq is ALREADY reverse-complemented onto the query strand and
    must not be complemented again; only the subject coordinate runs backwards.
    """
    sstep = 1 if sstart <= send else -1
    qpos, spos = qstart, sstart
    for qchar, schar in zip(qseq, sseq):
        q_gap, s_gap = qchar == "-", schar == "-"
        yield (None if q_gap else qpos), (None if s_gap else spos), qchar, schar
        if not q_gap:
            qpos += 1
        if not s_gap:
            spos += sstep


def base_at_query_pos(qseq: str, sseq: str, qstart: int, sstart: int, send: int, pos: int):
    """Read the subject base aligned to query position ``pos``.

    Walks the gapped alignment column by column instead of computing
    ``sstart + (pos - qstart)``. That arithmetic silently assumed an ungapped
    HSP, so a single indel anywhere between the HSP start and the diagnostic
    site shifted the read-out by one base -- the characteristic ONT error mode,
    and enough to flip a serotype call.

    Returns ``(base, contig_pos, strand)``. ``base`` is ``"-"`` when the
    subject carries a deletion at the site, or ``None`` if the alignment does
    not reach ``pos``.
    """
    strand = "+" if sstart <= send else "-"
    for qpos, spos, _qchar, schar in alignment_columns(qseq, sseq, qstart, sstart, send):
        if qpos == pos:
            return ("-" if schar == "-" else schar.upper()), spos, strand
    return None, None, strand


def triplet_at_query_pos(qseq: str, sseq: str, qstart: int, sstart: int, send: int, pos: int) -> dict:
    """Read the whole diagnostic codon -- query positions ``pos-2..pos``.

    Reading only the wobble base accepted any triplet ending in G as Trp: an
    AGG codon (Arg) in the subject was reported as TGG (Trp) and called as
    serotype 2 or 14. The codon is only interpretable if all three of its bases
    were actually recovered, contiguously, from one subject locus.

    Returns the read-out plus the subject coordinate of each base and a status
    naming exactly why it is or is not interpretable.
    """
    strand = "+" if sstart <= send else "-"
    sstep = 1 if strand == "+" else -1
    wanted = [pos - 2, pos - 1, pos]
    found: dict[int, tuple[str, int | None]] = {}

    for qpos, spos, _qchar, schar in alignment_columns(qseq, sseq, qstart, sstart, send):
        if qpos in wanted and qpos not in found:
            found[qpos] = ("-" if schar == "-" else schar.upper(), spos)
        if len(found) == 3:
            break

    chars = [found[p][0] if p in found else "?" for p in wanted]
    positions = [found[p][1] if p in found else None for p in wanted]
    triplet = "".join(chars)

    result = {"triplet": triplet, "positions": positions, "strand": strand,
              "contig_pos": positions[2], "base": chars[2],
              "covered": [p in found for p in wanted]}

    if "?" in chars:
        result["triplet_status"] = TRIPLET_INCOMPLETE
        return result
    if "-" in chars:
        result["triplet_status"] = TRIPLET_DELETED
        return result
    if any(c not in DNA for c in chars):
        result["triplet_status"] = TRIPLET_AMBIGUOUS
        return result
    # Contiguity: the three subject bases must be adjacent, in the strand's
    # direction. A gap in the query between them means the subject carries an
    # insertion inside the codon, so these three bases are not a codon of the
    # subject's own gene even though each aligns to the right query position.
    if any(positions[i + 1] - positions[i] != sstep for i in range(2)):
        result["triplet_status"] = TRIPLET_NONCONTIGUOUS
        return result
    result["triplet_status"] = TRIPLET_OK if triplet in RESOLVABLE_CODONS else TRIPLET_UNEXPECTED
    return result


def assess_coding_integrity(qseq: str, sseq: str, qstart: int, qend: int, qlen: int) -> dict:
    """Is the subject's copy of this reference CDS plausibly intact?

    Recovering a coordinate across an indel is not proof that the gene around
    it still reads. This looks for the two disruptions the alignment can
    actually show -- a frameshift, and a premature stop codon in the
    reference's reading frame -- and otherwise says plainly that it does not
    know.

      DISRUPTED   a frameshift or a premature stop was positively detected
      INTACT      the alignment spans the whole reference CDS and shows neither
      UNASSESSED  the alignment is partial, so the rest of the gene is unseen

    The reference is a complete CDS starting at query position 1, so query
    position p sits in codon ceil(p/3).
    """
    insertions = sum(1 for qc in qseq if qc == "-")
    deletions = sum(1 for sc in sseq if sc == "-")
    reasons = []
    if (insertions - deletions) % 3:
        reasons.append(f"frameshift:indel_balance={insertions - deletions}")

    # Project the subject onto query coordinates (insertion columns dropped),
    # then read complete codons of the reference frame.
    projected: dict[int, str] = {}
    qpos = qstart
    for qchar, schar in zip(qseq, sseq):
        if qchar == "-":
            continue
        projected[qpos] = "-" if schar == "-" else schar.upper()
        qpos += 1

    last_codon = qlen // 3  # the reference's own stop codon
    premature = []
    for codon_i in range(1, last_codon):
        p = codon_i * 3
        triplet = "".join(projected.get(x, "?") for x in (p - 2, p - 1, p))
        if triplet in STOP_CODONS:
            premature.append(codon_i)
    if premature:
        reasons.append(f"premature_stop:codon_{premature[0]}")

    spans_cds = min(qstart, qend) <= 1 and max(qstart, qend) >= qlen
    if reasons:
        status = CODING_DISRUPTED
    elif spans_cds:
        status = CODING_INTACT
    else:
        status = CODING_UNASSESSED
        reasons.append(f"alignment_covers_query_{min(qstart, qend)}-{max(qstart, qend)}_of_{qlen}")

    return {"coding_status": status, "coding_detail": ";".join(reasons),
            "insertions": insertions, "deletions": deletions}


# A codon that was fully recovered: all three bases present, contiguous and
# unambiguous. Only these can disagree with one another -- a locus whose codon
# could not be read says nothing about whether it agrees.
DETERMINED_TRIPLET_STATES = frozenset({TRIPLET_OK, TRIPLET_UNEXPECTED})


def codon_was_read(ev: dict) -> bool:
    return ev.get("triplet_status") in DETERMINED_TRIPLET_STATES


def collapse_to_loci(evidence: list[dict]) -> list[dict]:
    """Collapse alignments of several references to ONE physical subject locus.

    The four resolver references are 98-99.5% identical, so a single copy of
    cpsK in the assembly is hit by two or three of them. Different query names
    are therefore not evidence of multiple copies; overlapping subject
    coordinates on the same contig are one copy however many references found
    it. Each returned locus carries its best-scoring alignment plus the list of
    references that reached it.
    """
    by_contig = defaultdict(list)
    for ev in evidence:
        by_contig[ev["contig"]].append(ev)

    loci = []
    for contig, evs in by_contig.items():
        clusters: list[list[dict]] = []
        for ev in sorted(evs, key=lambda e: (e["s_lo"], e["s_hi"])):
            if clusters and ev["s_lo"] <= max(e["s_hi"] for e in clusters[-1]):
                clusters[-1].append(ev)
            else:
                clusters.append([ev])
        for cluster in clusters:
            best = max(cluster, key=lambda e: (e["bitscore"], e["ref_id"]))
            loci.append({**best,
                         "locus": f"{contig}:{min(e['s_lo'] for e in cluster)}-"
                                  f"{max(e['s_hi'] for e in cluster)}",
                         "ref_ids": sorted({e["ref_id"] for e in cluster}),
                         "n_alignments": len(cluster)})
    # A locus whose codon was actually read outranks a higher-scoring one whose
    # codon the alignment did not cover. On a fragmented assembly the longest
    # alignment is routinely the one that stops just short of the site, and
    # letting it become the primary turned a readable call into a no-call.
    # The trailing key makes the order total, so the reported coordinates do
    # not depend on the order BLAST happened to emit its rows in.
    return sorted(loci, key=lambda e: (not codon_was_read(e), -e["bitscore"],
                                       e["contig"], e["s_lo"]))


def implied_serotype(ev: dict) -> str | None:
    """The serotype one locus's evidence implies, or None if it implies none."""
    if ev.get("triplet_status") != TRIPLET_OK:
        return None
    meta = parse_resolver_meta(ev["ref_id"])
    triplet = ev["triplet"]
    if triplet in TRP_CODONS:
        return meta.get("G_serotype")
    if triplet in CYS_CODONS:
        return meta.get("CT_serotype")
    return None


def stage2_resolver_call(assembly_fa: str, resolver_refs_fa: str, threads: int, run_dir: Path,
                         config: dict, allowed_pair: str | None = None):
    """Read the diagnostic codon from every qualifying locus in the assembly.

    Returns None when nothing qualified. Otherwise returns the best-scoring
    locus's evidence, with ``loci`` carrying every distinct physical locus and
    ``conflict`` set when they do not agree. Keeping only the highest-bit-score
    HSP -- what this used to do -- silently picked a winner when two real copies
    of cpsK disagreed.
    """
    ensure_tool("blastn"); ensure_tool("makeblastdb")
    db_prefix = make_db_if_needed(assembly_fa, config["tmp_dir"])
    tsv_text = run_blast(resolver_refs_fa, db_prefix, threads, RESOLVER_OUTFMT)

    if config["keep_debug"]:
        stage2_tsv = run_dir / "resolver_vs_asm.tsv"
        stage2_tsv.write_text(tsv_text + ("\n" if tsv_text else ""))
        if config["gzip_debug"]:
            gzip_file(stage2_tsv)

    evidence = []
    for line in filter(None, tsv_text.splitlines()):
        parts = line.split("\t")
        if len(parts) < 13: continue
        (qseqid, sseqid, pident, length, qlen, evalue, bitscore,
         qstart, qend, sstart, send, qseq, sseq) = parts[:13]
        pident, length, qlen = float(pident), int(length), int(qlen)
        qstart, qend, sstart, send, bitscore = int(qstart), int(qend), int(sstart), int(send), float(bitscore)
        meta = parse_resolver_meta(qseqid)
        if allowed_pair and meta["pair"] != allowed_pair: continue
        pos = meta["pos"]
        q_lo, q_hi = min(qstart, qend), max(qstart, qend)
        # An HSP that does not reach the diagnostic base carries no information
        # about it at all. One that reaches it but not the rest of the codon is
        # kept and reported as INCOMPLETE rather than silently dropped.
        if not (q_lo <= pos <= q_hi): continue
        if pident < config["min_res_pid"] or length < config["min_res_alen"]: continue

        tri = triplet_at_query_pos(qseq, sseq, qstart, sstart, send, pos)
        coding = assess_coding_integrity(qseq, sseq, qstart, qend, qlen)
        evidence.append({"ref_id": qseqid, "contig": sseqid, "pair": meta["pair"],
                         "pident": pident, "length": length, "bitscore": bitscore,
                         "qstart": qstart, "qend": qend, "qlen": qlen,
                         "sstart": sstart, "send": send,
                         "s_lo": min(sstart, send), "s_hi": max(sstart, send),
                         **tri, **coding})

    if not evidence:
        return None

    loci = collapse_to_loci(evidence)
    for ev in loci:
        ev["implied_serotype"] = implied_serotype(ev)

    # Two loci conflict when both had their codon fully read and what the
    # codons MEAN differs: TGG against TGT, or TGG against AGG. Comparing the
    # codon spellings instead reported TGT against TGC as a conflict, although
    # both encode Cys161 and so both imply the same serotype. A codon outside
    # the scheme stands for itself, so it disagrees with any documented one.
    #
    # A locus whose codon could not be read (a partial alignment over a contig
    # boundary, an ambiguity code, a deletion) is recorded but does not
    # manufacture a disagreement: on a fragmented assembly the same physical
    # gene routinely produces one complete and one truncated alignment.
    readings = {ev["implied_serotype"] or ev["triplet"] for ev in loci if codon_was_read(ev)}
    conflict = len(readings) > 1

    primary = loci[0]
    return {**primary, "loci": loci, "n_loci": len(loci), "conflict": conflict,
            "conflicting_readings": sorted(readings) if conflict else []}


def resolver_status(ev: dict | None) -> str:
    """The Stage-2 outcome, as reported in the summary."""
    if ev is None:
        return S2_NO_HSP
    if ev.get("conflict"):
        return S2_CONFLICT
    if ev.get("triplet_status") != TRIPLET_OK:
        return S2_INVALID
    if ev.get("coding_status") == CODING_DISRUPTED:
        return S2_CODING_DISRUPTED
    return S2_OK


def interpret_resolver(ev: dict | None, config: dict) -> str | None:
    """The exact serotype the resolver evidence supports, or None.

    An exact call needs all of: no conflict between physical copies, a
    complete contiguous diagnostic codon, and that codon being one of the
    three the scheme defines. A positively-detected coding disruption also
    withholds it, unless `withhold_on_coding_disruption` is switched off.
    Anything short of that is withheld -- the family-level result is still
    reported by the caller.
    """
    if ev is None:
        return None
    if ev.get("conflict"):
        return None
    if ev.get("coding_status") == CODING_DISRUPTED \
            and (config or {}).get("withhold_on_coding_disruption", True):
        return None
    return implied_serotype(ev)


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
