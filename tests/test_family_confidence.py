"""Confidence BETWEEN families, established before resolving WITHIN one.

Stage 2 reads a single cpsK site that separates 2 from 1/2, or 1 from 14. It
answers no other question. So it may only be consulted once the evidence has
already placed the isolate in that family.

The tool used to reach Stage 2 whenever the top OR the runner-up individual
label happened to be a member of a resolvable family. With top="9" and
second="2" it selected family 2_vs_1_2, read the cpsK site, and overwrote the
unresolved "9" with a confident "2" -- a serotype the site was never asked
about.

Families are {1,14}, {2,1/2}, and every other type on its own. Within a
family, the score is the best score per marker class, not the sum: serotype 1
and serotype 14 references differ at a single site, so counting both would
manufacture two pieces of evidence out of one.
"""
import pytest

from helpers import REF_2_12, SUIS, main_config, resolver_double, run_stage1, stage1_double
from helpers import blast_row as row, stage1_cfg
from swineotype.stages import family_of


# --- family membership -------------------------------------------------

@pytest.mark.parametrize("serotype, family", [
    ("1", "1_vs_14"), ("14", "1_vs_14"),
    ("2", "2_vs_1_2"), ("1/2", "2_vs_1_2"),
    ("9", "type:9"), ("27", "type:27"), (None, None),
])
def test_family_of(serotype, family):
    assert family_of(serotype, stage1_cfg()) == family


# --- aggregation -------------------------------------------------------

def test_near_duplicate_references_in_one_family_are_not_independent_evidence():
    """serotype 1 and 14 differ at one site; both wzy references match one
    gene. Summing them would double the family's apparent support."""
    r = run_stage1({"wzy_1": ("1", "wzy"), "wzy_14": ("14", "wzy")},
               [row("wzy_1", 1800), row("wzy_14", 1810)])

    assert r["family_scores"]["1_vs_14"] == 1810.0, "best per marker class, not the sum"
    assert r["scores"]["1"] == 1800.0 and r["scores"]["14"] == 1810.0


def test_marker_classes_within_a_family_do_add():
    """wzx and wzy are different genes, so they are independent evidence."""
    r = run_stage1({"wzy_2": ("2", "wzy"), "wzx_2": ("2", "wzx")},
               [row("wzy_2", 1800), row("wzx_2", 1200)])
    assert r["family_scores"]["2_vs_1_2"] == 3000.0


# --- the five required decisions ---------------------------------------

def test_an_ordinary_within_family_tie_is_still_resolvable():
    """2 versus 1/2 at nearly equal scores is the normal case Stage 2 exists
    for. It must remain decisive at family level."""
    r = run_stage1({"wzy_2": ("2", "wzy"), "wzy_1_2": ("1/2", "wzy")},
               [row("wzy_2", 2000), row("wzy_1_2", 1990)])

    assert r["family_top"] == "2_vs_1_2"
    assert r["family_decisive"] is True
    assert r["family_fraction"] == pytest.approx(1.0)


def test_top_9_second_2_does_not_select_the_type_2_family():
    """Regression: the headline failure. A dominant serotype 9 with a weak
    serotype 2 runner-up used to hand Stage 2 the 2_vs_1_2 pair."""
    r = run_stage1({"wzy_9": ("9", "wzy"), "wzy_2": ("2", "wzy")},
               [row("wzy_9", 2000), row("wzy_2", 500)])

    assert r["top"] == "9" and r["second"] == "2"
    assert r["family_top"] == "type:9"
    assert r["family_decisive"] is True


def test_near_tied_competing_families_stay_unresolved():
    r = run_stage1({"wzy_14": ("14", "wzy"), "wzy_2": ("2", "wzy")},
               [row("wzy_14", 1000), row("wzy_2", 980)])

    assert set(r["family_scores"]) == {"1_vs_14", "2_vs_1_2"}
    assert r["family_delta"] == pytest.approx(20.0)
    assert r["family_decisive"] is False


def test_exact_markers_from_both_families_are_ambiguous():
    """A wzy from {1,14} and a wzy from {2,1/2}, both strong. The cpsK site
    cannot arbitrate between families, so nothing may be resolved."""
    r = run_stage1({"wzy_14": ("14", "wzy"), "wzy_2": ("2", "wzy")},
               [row("wzy_14", 2000), row("wzy_2", 2000)])

    assert r["family_delta"] == 0.0
    assert r["family_fraction"] == pytest.approx(0.5)
    assert r["family_decisive"] is False


def test_a_competitor_is_not_hidden_by_two_same_family_front_runners():
    """The first two INDIVIDUAL labels are 2 and 1/2 -- one family -- so the
    old top-two check saw no disagreement at all. Serotype 9, only 10 bits
    behind the family, was never considered."""
    r = run_stage1({"wzy_2": ("2", "wzy"), "wzy_1_2": ("1/2", "wzy"), "wzy_9": ("9", "wzy")},
               [row("wzy_2", 1000), row("wzy_1_2", 995), row("wzy_9", 990)])

    assert [r["top"], r["second"]] == ["2", "1/2"], "both front-runners are one family"
    assert r["family_top"] == "2_vs_1_2"
    assert r["family_second"] == "type:9", "the real competitor is visible"
    assert r["family_delta"] == pytest.approx(10.0)
    assert r["family_decisive"] is False


# --- what process_one does with it -------------------------------------

def test_unresolved_families_produce_no_call_and_no_stage_2(patched_stages, tmp_path):
    s1 = stage1_double(top="9", second="2", family_top="type:9",
                       family_second="2_vs_1_2", family_decisive=False,
                       family_fraction=0.52, family_delta=40.0)
    main_mod = patched_stages(s1, forbid_stage2=True)

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "NO_CALL_FAMILY_AMBIGUOUS"
    assert r["final_serotype"] == ""
    assert r["family_serotype"] == ""
    assert "competing_cps_families:type:9/2_vs_1_2" in r["warnings"]


def test_a_decisive_singleton_family_is_called_without_stage_2(patched_stages, tmp_path):
    s1 = stage1_double(top="9", second="2", family_top="type:9",
                       family_second="2_vs_1_2")
    main_mod = patched_stages(s1, forbid_stage2=True)

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "STAGE1"
    assert r["final_serotype"] == "9", "not 2 -- the cpsK site is never consulted"


def test_a_decisive_resolvable_family_reaches_stage_2(patched_stages, tmp_path):
    s1 = stage1_double(top="2", second="1/2", family_top="2_vs_1_2")
    main_mod = patched_stages(s1, s2=resolver_double(REF_2_12, "TGG"))

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "STAGE2"
    assert r["final_serotype"] == "2"


def test_family_result_survives_a_failed_exact_resolution(patched_stages, tmp_path):
    """Withholding the exact member must not throw away the family."""
    s1 = stage1_double(top="2", second="1/2", family_top="2_vs_1_2")
    main_mod = patched_stages(s1, s2=resolver_double(REF_2_12, "AGG"))

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "FAMILY_ONLY"
    assert r["final_serotype"] == ""
    assert r["family_serotype"] == "2 or 1/2"
    assert r["stage1_family"] == "2_vs_1_2"


def test_a_stage1_result_with_no_family_assessment_is_not_assumed_confident(patched_stages, tmp_path):
    """Defensive: an old-shaped Stage-1 result has not been assessed between
    families, so it must not be treated as though it had been."""
    legacy = {"top": "2", "second": "1/2", "decisive": True, "delta": 5000.0,
              "must_stage2_for_pair": True, "top_species": SUIS, "top_wzx_only": None}
    main_mod = patched_stages(legacy, forbid_stage2=True)

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))
    assert r["status"] == "NO_CALL_FAMILY_AMBIGUOUS"


def test_the_family_decision_picks_the_resolver_pair_not_the_top_label():
    """A family whose members each contribute a different marker class can
    out-score the single highest individual type. The family that was actually
    assessed for confidence is the one Stage 2 must be pointed at."""
    r = run_stage1({"wzy_9": ("9", "wzy"),
                "wzy_2": ("2", "wzy"), "wzx_1_2": ("1/2", "wzx")},
               [row("wzy_9", 1000), row("wzy_2", 900), row("wzx_1_2", 900)])

    assert r["top"] == "9", "the highest individual type"
    assert r["family_top"] == "2_vs_1_2", "but not the best-supported family"
    assert r["family_scores"]["2_vs_1_2"] == 1800.0


def test_a_family_top_that_disagrees_with_the_top_label_is_recorded(patched_stages, tmp_path):
    s1 = stage1_double(top="9", second="2", family_top="2_vs_1_2",
                       family_second="type:9")
    main_mod = patched_stages(s1, s2=resolver_double(REF_2_12, "TGG"))

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "STAGE2"
    assert r["final_serotype"] == "2"
    assert "family_disagrees_with_top_label:2_vs_1_2/None" not in r["warnings"]
    assert "family_disagrees_with_top_label" not in r["warnings"], \
        "top='9' implies no pair at all, so there is nothing to disagree with"


def test_a_top_label_pointing_at_the_other_family_is_flagged(patched_stages, tmp_path):
    s1 = stage1_double(top="14", second="2", family_top="2_vs_1_2",
                       family_second="1_vs_14")
    main_mod = patched_stages(s1, s2=resolver_double(REF_2_12, "TGG"))

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert "family_disagrees_with_top_label:2_vs_1_2/1_vs_14" in r["warnings"]
    assert r["final_serotype"] == "2", "the family decision governs"
