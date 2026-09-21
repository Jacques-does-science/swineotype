"""Shared constants and builders for the swineotype test suite.

Imported as a top-level module: pytest prepends the tests/ directory to
sys.path, and `tests` itself is a name already taken by an unrelated package
in some environments.
"""
from __future__ import annotations

SUIS = "Streptococcus suis"

# Real reference conventions, so nobody reads these fixtures and infers the
# wrong polarity: TGG (Trp161) -> Gal -> serotype 2 or 14.
REF_1_14 = "cps14K|pair=1_vs_14|pos=492|G_serotype=14|CT_serotype=1"
REF_1_14_ALT = "cps1L|pair=1_vs_14|pos=492|G_serotype=14|CT_serotype=1"
REF_2_12 = "cps2K|pair=2_vs_1_2|pos=483|G_serotype=2|CT_serotype=1/2"
REF_2_12_ALT = "cps1/2K|pair=2_vs_1_2|pos=483|G_serotype=2|CT_serotype=1/2"

RESOLVABLE = ("1_vs_14", "2_vs_1_2")


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
    s1 = {"top": top, "second": second, "decisive": True, "delta": 5000.0,
          "fraction": 1.0, "must_stage2_for_pair": family_top in RESOLVABLE,
          "top_species": SUIS, "top_wzx_only": None,
          "family_top": family_top, "family_second": family_second,
          "family_decisive": family_decisive,
          "family_fraction": 1.0, "family_delta": 5000.0,
          "family_types": {}, "family_scores": {}}
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
          "conflicting_serotypes": []}
    ev.update(over)
    from swineotype.stages import implied_serotype
    ev.setdefault("implied_serotype", implied_serotype(ev))
    # A copy, not `ev` itself: a self-referential dict cannot be serialised.
    ev.setdefault("loci", [dict(ev)])
    return ev
