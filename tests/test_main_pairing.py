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
    # Fall back to the runner-up only when the top hit is in neither pair.
    ("9", "2", "2_vs_1_2"),
    ("9", "14", "1_vs_14"),
    ("9", "7", None),
    (None, None, None),
])
def test_choose_pair_prefers_the_top_hit(top, second, expected):
    assert choose_pair(top, second, CONFIG) == expected


def test_choose_pair_never_returns_the_runner_ups_pair_over_the_tops():
    """Exhaustive: whenever the top hit belongs to a pair, that pair wins."""
    everything = ["1", "14", "2", "1/2", "9", None]
    for top in everything:
        for second in everything:
            result = choose_pair(top, second, CONFIG)
            if pair_of(top, CONFIG):
                assert result == pair_of(top, CONFIG), f"top={top} second={second}"
