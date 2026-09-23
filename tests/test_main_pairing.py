"""Regression tests for resolver-pair selection.

Stage 1 places an isolate in the {1,14} or the {2,1/2} locus group; Stage 2
then reads the cpsK site using only that group's references. Because the four
references are 98-99.5% identical to one another, the pair filter is the ONLY
thing preventing a confident cross-pair call -- percent-identity thresholds
cannot separate them. So pair selection has to be exactly right.

These used to exercise choose_pair(top, second). That function ignored its
second argument once the runner-up fallback was removed, so it was deleted;
the cases now run through process_one(), where a wrong pair actually does harm.
"""
import pytest

from helpers import main_config, stage1_double
from swineotype.stages import pair_of

CONFIG = {"pair_1_14": {"1", "14"}, "pair_2_1_2": {"2", "1/2"}}


@pytest.mark.parametrize("serotype, expected", [
    ("1", "1_vs_14"), ("14", "1_vs_14"),
    ("2", "2_vs_1_2"), ("1/2", "2_vs_1_2"),
    ("9", None), (None, None),
])
def test_pair_of(serotype, expected):
    assert pair_of(serotype, CONFIG) == expected


@pytest.fixture
def stage2_pair(monkeypatch, tmp_path):
    """Run process_one for a (top, second) Stage-1 result; return the pair
    Stage 2 was asked to read, or None if Stage 2 never ran."""
    from swineotype import main as main_mod

    def run(top, second):
        asked = {}

        def stage2(assembly, refs, threads, run_dir, config, allowed_pair):
            asked["pair"] = allowed_pair

        monkeypatch.setattr(main_mod, "stage1_score",
                            lambda *a, **k: stage1_double(top=top, second=second))
        monkeypatch.setattr(main_mod, "ensure_unix_line_endings", lambda p, t: p)
        monkeypatch.setattr(main_mod, "stage2_resolver_call", stage2)
        row = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))
        return asked.get("pair"), row

    return run


@pytest.mark.parametrize("top, second, expected", [
    # The top hit decides, whichever pair the runner-up belongs to.
    ("2", "1/2", "2_vs_1_2"),
    ("1", "14", "1_vs_14"),
    # Regression: top in {2,1/2} with a runner-up in {1,14}. The old code
    # tested pair_1_14 against top AND second before ever testing pair_2_1_2,
    # so this returned "1_vs_14" and a true serotype 2 was reported as 14.
    ("2", "14", "2_vs_1_2"),
    ("2", "1", "2_vs_1_2"),
    ("1/2", "1", "2_vs_1_2"),
    # Mirror case, which the old code happened to get right.
    ("14", "2", "1_vs_14"),
    ("1", "1/2", "1_vs_14"),
])
def test_stage_2_reads_the_top_hits_pair(stage2_pair, top, second, expected):
    pair, _ = stage2_pair(top, second)
    assert pair == expected


@pytest.mark.parametrize("top, second", [
    # Regression, and the expectation this file used to encode backwards: a
    # top hit outside both resolvable families selects NO pair. It used to
    # fall back to the runner-up, so top="9" second="2" chose 2_vs_1_2, Stage
    # 2 read the cpsK site and overwrote the unresolved 9 with a confident
    # "2". The cpsK site distinguishes 2 from 1/2; it is silent on whether the
    # isolate is a 9, so reading it here invents confidence.
    ("9", "2"), ("9", "14"), ("9", "1/2"), ("9", "1"), ("9", "7"),
])
def test_a_top_hit_outside_both_pairs_never_reaches_stage_2(stage2_pair, top, second):
    pair, row = stage2_pair(top, second)
    assert pair is None, "the runner-up's pair must not be consulted"
    assert row["final_serotype"] == top
    assert row["status"] == "STAGE1"
