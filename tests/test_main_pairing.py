"""Regression tests for resolver-pair selection.

Stage 1 places an isolate in the {1,14} or the {2,1/2} locus group; Stage 2
then reads the cpsK site using only that group's references. Because the four
references are 98-99.5% identical to one another, the pair filter is the ONLY
thing preventing a confident cross-pair call -- percent-identity thresholds
cannot separate them. So pair selection has to be exactly right.
"""
import pytest

from swineotype.main import choose_pair, pair_of

CONFIG = {"pair_1_14": {"1", "14"}, "pair_2_1_2": {"2", "1/2"}}


@pytest.mark.parametrize("serotype, expected", [
    ("1", "1_vs_14"), ("14", "1_vs_14"),
    ("2", "2_vs_1_2"), ("1/2", "2_vs_1_2"),
    ("9", None), (None, None),
])
def test_pair_of(serotype, expected):
    assert pair_of(serotype, CONFIG) == expected


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
    # Regression, and the expectation this file used to encode backwards:
    # a top hit outside both resolvable families selects NO pair. It used to
    # fall back to the runner-up, so top="9" second="2" chose 2_vs_1_2, Stage
    # 2 read the cpsK site and overwrote the unresolved 9 with a confident
    # "2". The cpsK site distinguishes 2 from 1/2; it is silent on whether the
    # isolate is a 9, so reading it here invents confidence.
    ("9", "2", None),
    ("9", "14", None),
    ("9", "7", None),
    (None, "2", None),
    (None, None, None),
])
def test_choose_pair_uses_only_the_top_hit(top, second, expected):
    assert choose_pair(top, second, CONFIG) == expected


def test_choose_pair_never_returns_the_runner_ups_pair_over_the_tops():
    """Exhaustive: the runner-up never influences the answer."""
    everything = ["1", "14", "2", "1/2", "9", None]
    for top in everything:
        for second in everything:
            assert choose_pair(top, second, CONFIG) == pair_of(top, CONFIG), \
                f"top={top} second={second}"
