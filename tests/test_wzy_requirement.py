"""A serotype call requires its serotype-specific gene (wzy).

wzx (flippase) is conserved across serotypes; a single wzx routinely matches
four different references at 95-98% identity and full coverage. wzy
(polymerase) is the serotype-specific gene -- which is why the published
multiplex PCR schemes target it.

Scoring them equally meant a conserved wzx alone could crown a serotype. A real
isolate carrying a novel capsular locus reported `stage1_top=27` on the
strength of its wzx, while its own wzy matched nothing in the panel at usable
identity (best: 66.8% to serotype 16). The right answer there is "this locus is
not in the panel", not a serotype.

Requiring wzx as well would be wrong: serotype 14 has no wzx reference.
"""
import pytest
from unittest.mock import patch
from pathlib import Path

from helpers import SUIS, main_config, stage1_double
from swineotype.config import load_config
from swineotype.stages import parse_whitelist_headers, stage1_score


def cfg(**over):
    c = {"min_pid": 85.0, "min_cov": 0.8, "plurality": 0.6, "delta": 100,
         "ambig_set": {"1", "14", "2", "1/2"}, "keep_debug": False,
         "tmp_dir": "tmp", "require_wzy": 1,
         "pair_1_14": {"1", "14"}, "pair_2_1_2": {"2", "1/2"}}
    c.update(over)
    return c


def row(allele, bits, cov_len=100, qlen=100, pid=100, contig="c1"):
    """One stage-1 BLAST row (11 columns)."""
    return "\t".join([allele, contig, str(pid), str(cov_len), str(qlen), "0",
                      str(bits), "1", str(cov_len), "1", str(cov_len)])


def run_stage1(alleles, blast_rows, config=None, species=None):
    """alleles: {allele_id: (type, geneclass)}"""
    a2t = {a: t for a, (t, g) in alleles.items()}
    a2g = {a: g for a, (t, g) in alleles.items()}
    t2s = species or {t: SUIS for t, _ in alleles.values()}
    with patch("swineotype.stages.ensure_tool"), \
         patch("swineotype.stages.make_db_if_needed", return_value="db"), \
         patch("swineotype.stages.run_blast", return_value="\n".join(blast_rows)), \
         patch("swineotype.stages.parse_whitelist_headers", return_value=(a2t, a2g, t2s)):
        return stage1_score("a.fasta", "w.fasta", 1, Path("run"), config or cfg())


# --- data invariant ----------------------------------------------------

def test_every_type_in_the_panel_has_a_wzy_reference():
    """require_wzy can only be safe if no type depends on wzx alone."""
    a2t, a2g, _ = parse_whitelist_headers(load_config()["wzxwzy_fasta"])
    by_type = {}
    for allele, t in a2t.items():
        by_type.setdefault(t, set()).add(a2g.get(allele))
    missing = sorted(t for t, genes in by_type.items() if "wzy" not in genes)
    assert not missing, f"types with no wzy reference would become uncallable: {missing}"


def test_serotype_14_has_no_wzx_reference():
    """Documents why the rule is wzy-only rather than 'both genes'."""
    a2t, a2g, _ = parse_whitelist_headers(load_config()["wzxwzy_fasta"])
    genes_14 = {a2g.get(a) for a, t in a2t.items() if t == "14"}
    assert genes_14 == {"wzy"}


# --- the rule ----------------------------------------------------------

def test_wzx_only_type_cannot_be_called():
    """Regression: the real-isolate failure.

    A conserved wzx matches serotype 27 best, but no wzy passes. 27 must not
    be reported as the call; it is reported as the nearest wzx relative.
    """
    r = run_stage1(
        {"wzx_27": ("27", "wzx"), "wzy_27": ("27", "wzy")},
        [row("wzx_27", 2407)],                      # wzy_27 absent = rejected
    )
    assert r["top"] is None, "a wzx-only match must not produce a call"
    assert r["top_wzx_only"] == "27", "but it should be reported as a lead"
    assert r["scores"]["27"] == 2407.0, "raw score is still recorded for debug"


def test_wzy_support_makes_a_type_callable():
    r = run_stage1(
        {"wzx_9": ("9", "wzx"), "wzy_9": ("9", "wzy")},
        [row("wzx_9", 1200), row("wzy_9", 1300)],
    )
    assert r["top"] == "9"
    assert r["top_wzx_only"] is None
    assert r["genes_by_type"]["9"] == ["wzx", "wzy"]


def test_serotype_14_is_callable_from_wzy_alone():
    """The whole reason the rule is wzy-only and not 'both'."""
    r = run_stage1({"wzy_14": ("14", "wzy")}, [row("wzy_14", 2000)])
    assert r["top"] == "14"


def test_wzy_supported_type_beats_a_higher_scoring_wzx_only_type():
    """This is the exact shape of the failing isolate: the wzx-only type
    outscores the real one, and must still lose."""
    r = run_stage1(
        {"wzx_27": ("27", "wzx"), "wzy_27": ("27", "wzy"),
         "wzx_9": ("9", "wzx"), "wzy_9": ("9", "wzy")},
        [row("wzx_27", 2407), row("wzx_9", 900), row("wzy_9", 800)],
    )
    assert r["top"] == "9"
    assert r["top_wzx_only"] == "27"


# --- plurality is measured among callable candidates -------------------

def test_plurality_ignores_wzx_only_cross_hits():
    """Regression: serotype 27 used to return NO_CALL.

    Its wzx cross-hybridises with serotypes 1, 2 and 1/2, so those three
    landed in the plurality denominator and dragged the fraction to 0.41
    against a 0.60 gate -- despite a delta of 2368 against a threshold of 100.
    """
    alleles = {"wzx_27": ("27", "wzx"), "wzy_27": ("27", "wzy"),
               "wzx_1": ("1", "wzx"), "wzy_1": ("1", "wzy"),
               "wzx_2": ("2", "wzx"), "wzy_2": ("2", "wzy")}
    rows = [row("wzx_27", 2300), row("wzy_27", 2300),   # the real serotype
            row("wzx_1", 2269), row("wzx_2", 2273)]     # wzx cross-hits only

    strict = run_stage1(alleles, rows, cfg(require_wzy=1))
    assert strict["top"] == "27"
    assert strict["fraction"] == pytest.approx(1.0)
    assert strict["decisive"] is True

    legacy = run_stage1(alleles, rows, cfg(require_wzy=0))
    assert legacy["fraction"] < 0.6, "old behaviour: diluted by wzx-only types"
    assert legacy["decisive"] is False


def test_require_wzy_can_be_switched_off():
    r = run_stage1(
        {"wzx_27": ("27", "wzx"), "wzy_27": ("27", "wzy")},
        [row("wzx_27", 2407)],
        cfg(require_wzy=0),
    )
    assert r["top"] == "27", "with the rule disabled, old behaviour returns"


# --- the reported result ------------------------------------------------

def test_no_wzy_match_reports_a_lead_without_claiming_an_intact_locus(patched_stages, tmp_path):
    """A wzx-only match establishes that no qualifying reference wzy was
    found. It does NOT establish an intact locus, and it does not establish
    novelty -- an absent, truncated or contig-broken wzy looks the same. The
    old `candidate_novel_capsular_locus` warning asserted both."""
    s1 = stage1_double(top=None, family_top=None, family_decisive=False,
                       decisive=False, delta=0.0, top_species=None, top_wzx_only="27")
    main_mod = patched_stages(s1, forbid_stage2=True)

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "NO_WZY_MATCH"
    assert r["final_serotype"] == ""
    assert r["stage1_top"] == "", "no serotype may be implied"
    assert r["stage1_wzx_only"] == "27", "the lead is reported separately"
    assert r["matched_reference_taxon"] == "", "wzx is not species-discriminating either"
    assert r["species"] == ""
    assert "no_reference_wzy_matched" in r["warnings"]
    assert "nearest_wzx_relative:cps_type_27" in r["warnings"]
    assert "locus_integrity_not_assessed" in r["warnings"]
    assert "novel" not in r["warnings"], "novelty was never established"


def test_no_cps_match_infers_no_species_at_all(patched_stages, tmp_path):
    """Regression: `species or config["target_species"]` filled the column
    from a default, so an assembly that matched nothing -- an unrelated
    organism, an empty file -- was reported as Streptococcus suis."""
    s1 = stage1_double(top=None, family_top=None, family_decisive=False,
                       decisive=False, delta=0.0, top_species=None, top_wzx_only=None)
    main_mod = patched_stages(s1, forbid_stage2=True)

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "NO_CPS_MATCH"
    assert r["species"] == ""
    assert r["matched_reference_taxon"] == ""
    assert r["input_species"] == ""
    assert r["species_assessment"] == "NOT_ASSESSED"
    assert r["final_serotype"] == ""
    assert r["family_serotype"] == ""
    assert "no_cps_reference_matched" in r["warnings"]


def test_stage2_status_is_surfaced(patched_stages, tmp_path):
    """It was computed in process_one and then discarded, so a user could not
    see why Stage 2 produced nothing."""
    s1 = stage1_double(top="2", second="1/2", family_top="2_vs_1_2",
                       family_second=None, decisive=False, delta=10.0)
    main_mod = patched_stages(s1, s2=None)

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))
    assert r["stage2_status"] == "NO_HSP_OR_LOW_QUAL"
    assert r["status"] == "FAMILY_ONLY"
    assert r["final_serotype"] == ""
    assert r["family_serotype"] == "2 or 1/2", "the family result survives"


def test_stage2_status_ok_on_a_successful_resolve(patched_stages, tmp_path):
    from helpers import REF_2_12, resolver_double
    s1 = stage1_double(top="2", second="1/2", family_top="2_vs_1_2",
                       decisive=False, delta=10.0)
    main_mod = patched_stages(s1, s2=resolver_double(REF_2_12, "TGG", contig_pos=883))

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))
    assert r["stage2_status"] == "OK"
    assert r["final_serotype"] == "2"
    assert r["triplet"] == "TGG"
