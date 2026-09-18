"""Invariants over the real data/ files -- no mocks.

The cpsK polarity was inverted across three commits (#14 -> #15 -> #16) and
ended up correct with nothing locking it in; the one test that touched it was
red and asserted the opposite convention. These tests assert the published
biology directly against the shipped references.

Reference: a single amino-acid polymorphism at CpsK residue 161 sets whether
the glycosyltransferase adds Gal or GalNAc to the CPS side chain.
  Trp161 (codon TGG) -> Gal    -> serotype 2  (or 14 in the 1/14 locus)
  Cys161 (codon TGY) -> GalNAc -> serotype 1/2 (or 1 in the 1/14 locus)
Roy et al. 2017, Sci Rep 7:4066, doi:10.1038/s41598-017-04403-3
Athey et al. 2016, BMC Microbiol 16:162, doi:10.1186/s12866-016-0782-8
"""
import pytest

from swineotype.config import load_config
from swineotype.stages import interpret_resolver, parse_resolver_meta

TRP = {"TGG"}
CYS = {"TGT", "TGC"}

# Expected call for each shipped reference at its own declared position.
EXPECTED = {
    "cps2K": ("2_vs_1_2", 483, "G", "2"),
    "cps1/2K": ("2_vs_1_2", 483, "T", "1/2"),
    "cps14K": ("1_vs_14", 492, "G", "14"),
    "cps1L": ("1_vs_14", 492, "T", "1"),
}


def read_fasta(path):
    seqs, header = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
                seqs[header] = []
            elif header:
                seqs[header].append(line)
    return {h: "".join(parts) for h, parts in seqs.items()}


@pytest.fixture(scope="module")
def resolver_refs():
    return read_fasta(load_config()["resolver_refs_fasta"])


def test_all_four_references_are_present(resolver_refs):
    names = {h.split("|")[0] for h in resolver_refs}
    assert names == set(EXPECTED)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_reference_header_matches_expected_pair_and_position(resolver_refs, name):
    header = next(h for h in resolver_refs if h.split("|")[0] == name)
    meta = parse_resolver_meta(header)
    exp_pair, exp_pos, _, _ = EXPECTED[name]
    assert meta["pair"] == exp_pair
    assert meta["pos"] == exp_pos


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_diagnostic_site_holds_the_expected_base(resolver_refs, name):
    header = next(h for h in resolver_refs if h.split("|")[0] == name)
    seq = resolver_refs[header]
    _, pos, exp_base, _ = EXPECTED[name]
    assert seq[pos - 1].upper() == exp_base


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_diagnostic_site_is_the_third_base_of_a_trp_or_cys_codon(resolver_refs, name):
    """The declared position must be the wobble base of codon 161.

    This is what makes the differing positions (483 vs 492) legitimate rather
    than a typo: the 1/14 references carry 9 extra bases at the 5' end.
    """
    header = next(h for h in resolver_refs if h.split("|")[0] == name)
    seq = resolver_refs[header].upper()
    _, pos, _, expected_serotype = EXPECTED[name]
    codon = seq[pos - 3:pos]
    assert codon in TRP | CYS, f"{name}: codon at {pos} is {codon!r}, not Trp/Cys"
    assert pos % 3 == 0, f"{name}: position {pos} is not a codon-third base"
    if codon in TRP:
        assert expected_serotype in ("2", "14"), "Trp161 -> Gal -> serotype 2 / 14"
    else:
        assert expected_serotype in ("1/2", "1"), "Cys161 -> GalNAc -> serotype 1/2 / 1"


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_interpret_resolver_returns_the_published_serotype(resolver_refs, name):
    header = next(h for h in resolver_refs if h.split("|")[0] == name)
    _, _, base, expected_serotype = EXPECTED[name]
    assert interpret_resolver({"ref_id": header, "base": base}, {}) == expected_serotype


def test_the_two_pairs_share_flanking_context(resolver_refs):
    """Sanity-check the 9 bp offset: the sequence around the diagnostic base is
    the same in both pairs, which is why the position differs but the residue
    does not."""
    contexts = set()
    for name, (_, pos, _, _) in EXPECTED.items():
        header = next(h for h in resolver_refs if h.split("|")[0] == name)
        seq = resolver_refs[header].upper()
        contexts.add(seq[pos - 9:pos - 1] + seq[pos:pos + 8])
    assert len(contexts) == 1, f"flanking context diverges between references: {contexts}"
