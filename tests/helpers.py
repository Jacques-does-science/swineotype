"""Shared constants and builders for the swineotype test suite.

Imported as a top-level module: pytest prepends the tests/ directory to
sys.path, and `tests` itself is a name already taken by an unrelated package
in some environments.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

SUIS = "Streptococcus suis"

# Real reference conventions, so nobody reads these fixtures and infers the
# wrong polarity: TGG (Trp161) -> Gal -> serotype 2 or 14.
REF_1_14 = "cps14K|pair=1_vs_14|pos=492|G_serotype=14|CT_serotype=1"
REF_1_14_ALT = "cps1L|pair=1_vs_14|pos=492|G_serotype=14|CT_serotype=1"
REF_2_12 = "cps2K|pair=2_vs_1_2|pos=483|G_serotype=2|CT_serotype=1/2"
REF_2_12_ALT = "cps1/2K|pair=2_vs_1_2|pos=483|G_serotype=2|CT_serotype=1/2"


def main_config(**over):
    """A config of the shape process_one() expects."""
    cfg = {"target_species": SUIS, "pair_1_14": {"1", "14"}, "pair_2_1_2": {"2", "1/2"},
           "require_wzy": 1, "wzxwzy_fasta": "w.fasta", "resolver_refs_fasta": "r.fasta"}
    cfg.update(over)
    return cfg


def stage1_double(top=None, second=None, family_top=None, family_second=None,
                  family_decisive=True, config=None, **over):
    """A Stage-1 result of the shape stage1_score() now returns.

    `family_decisive` defaults to True so a test that is not about family
    confidence need not think about it. process_one treats a result carrying
    no family fields at all as NOT assessed, which is the honest default for
    real data but noise in a focused unit test.
    """
    from swineotype.stages import family_of
    cfg = config or main_config()
    if family_top is None and top is not None:
        family_top = family_of(top, cfg)
    s1 = {"top": top, "second": second,
          "top_species": SUIS, "top_wzx_only": None,
          "family_top": family_top, "family_second": family_second,
          "family_decisive": family_decisive,
          "family_fraction": 1.0, "family_delta": 5000.0, "family_scores": {}}
    s1.update(over)
    return s1


def resolver_double(ref_id, triplet, *, contig="c1", contig_pos=900, strand="+",
                    triplet_status=None, coding_status="UNASSESSED", **over):
    """One Stage-2 locus's evidence, as stage2_resolver_call() returns it."""
    if triplet_status is None:
        triplet_status = "OK" if triplet in ("TGG", "TGT", "TGC") else "UNEXPECTED_CODON"
    ev = {"ref_id": ref_id, "contig": contig, "contig_pos": contig_pos,
          "strand": strand, "base": triplet[-1], "triplet": triplet,
          "triplet_status": triplet_status, "coding_status": coding_status,
          "coding_detail": "", "bitscore": 1800.0,
          "locus": f"{contig}:1-1000", "n_loci": 1, "conflict": False,
          "conflicting_readings": []}
    ev.update(over)
    from swineotype.stages import implied_serotype
    ev.setdefault("implied_serotype", implied_serotype(ev))
    # A copy, not `ev` itself: a self-referential dict cannot be serialised.
    ev.setdefault("loci", [dict(ev)])
    return ev


# --- a mocked stage1_score() ------------------------------------------

def stage1_cfg(**over):
    c = {"min_pid": 85.0, "min_cov": 0.8, "plurality": 0.6, "delta": 100,
         "keep_debug": False, "tmp_dir": "tmp", "require_wzy": 1,
         "pair_1_14": {"1", "14"}, "pair_2_1_2": {"2", "1/2"}}
    c.update(over)
    return c


def blast_row(allele, bits, length=100, qlen=100, pid=100, contig="c1"):
    """One stage-1 BLAST row (11 columns), aligned from query position 1."""
    return "\t".join([allele, contig, str(pid), str(length), str(qlen), "0",
                      str(bits), "1", str(length), "1", str(length)])


def run_stage1(alleles, blast_rows, config=None, species=None):
    """stage1_score() over canned BLAST rows. alleles: {allele_id: (type, geneclass)}"""
    from swineotype.stages import stage1_score
    a2t = {a: t for a, (t, g) in alleles.items()}
    a2g = {a: g for a, (t, g) in alleles.items()}
    t2s = species or {t: SUIS for t, _ in alleles.values()}
    with patch("swineotype.stages.make_db_if_needed", return_value="db"), \
         patch("swineotype.stages.run_blast", return_value="\n".join(blast_rows)), \
         patch("swineotype.stages.parse_whitelist_headers", return_value=(a2t, a2g, t2s)):
        return stage1_score("a.fasta", "w.fasta", 1, Path("run"), config or stage1_cfg())


# --- a mocked stage2_resolver_call() ----------------------------------

STAGE2_CFG = {"min_res_pid": 90, "min_res_alen": 100, "keep_debug": False, "tmp_dir": "tmp"}


def resolver_row(qseqid, pos, triplet, contig="c1", sstart=1, bitscore=2000,
                 qlen=1000, pident=100):
    """A 13-column resolver row whose subject carries `triplet` at pos-2..pos.

    The query carries the canonical TGG so the row is a plausible alignment;
    only the subject varies.
    """
    qseq = "A" * (pos - 3) + "TGG" + "A" * (qlen - pos)
    sseq = "A" * (pos - 3) + triplet + "A" * (qlen - pos)
    return "\t".join([qseqid, contig, str(pident), str(qlen), str(qlen), "0", str(bitscore),
                      "1", str(qlen), str(sstart), str(sstart + qlen - 1), qseq, sseq])


def run_stage2(rows, allowed_pair="2_vs_1_2"):
    """stage2_resolver_call() over canned resolver BLAST rows."""
    from swineotype.stages import stage2_resolver_call
    with patch("swineotype.stages.make_db_if_needed", return_value="db"), \
         patch("swineotype.stages.run_blast", return_value="\n".join(rows)):
        return stage2_resolver_call("a.fa", "r.fa", 1, Path("run"), STAGE2_CFG, allowed_pair)
