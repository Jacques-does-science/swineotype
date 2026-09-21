"""Invariants over the real data/ files -- no mocks, no network.

The cpsK polarity was inverted across three commits (#14 -> #15 -> #16) and
ended up correct with nothing locking it in; the one test that touched it was
red and asserted the opposite convention. Separately, two of the four shipped
sequences carried a two-base deletion relative to their source records, which
no test could see because nothing checked the sequence as a whole. These tests
assert the published biology AND the byte-level identity of the shipped
references against the pinned manifest.

Reference: a single amino-acid polymorphism at CpsK residue 161 sets whether
the glycosyltransferase adds Gal or GalNAc to the CPS side chain.
  Trp161 (codon TGG) -> Gal    -> serotype 2  (or 14 in the 1/14 locus)
  Cys161 (codon TGY) -> GalNAc -> serotype 1/2 (or 1 in the 1/14 locus)
Roy et al. 2017, Sci Rep 7:4066, doi:10.1038/s41598-017-04403-3
Athey et al. 2016, BMC Microbiol 16:162, doi:10.1186/s12866-016-0782-8
"""
import hashlib
import json
import sys
from pathlib import Path

import pytest

from swineotype.config import load_config
from swineotype.stages import interpret_resolver, parse_resolver_meta

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import resolver_refs  # noqa: E402

TRP = {"TGG"}
CYS = {"TGT", "TGC"}
STOPS = {"TAA", "TAG", "TGA"}

# Expected call for each shipped reference at its own declared position.
EXPECTED = {
    "cps2K": ("2_vs_1_2", 483, "TGG", "2"),
    "cps1/2K": ("2_vs_1_2", 483, "TGT", "1/2"),
    "cps14K": ("1_vs_14", 492, "TGG", "14"),
    "cps1L": ("1_vs_14", 492, "TGT", "1"),
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
def resolver_refs_fasta():
    return read_fasta(load_config()["resolver_refs_fasta"])


@pytest.fixture(scope="module")
def manifest():
    with open(Path(load_config()["data_dir"]) / "suis_resolver_refs.manifest.json") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def manifest_by_id(manifest):
    return {e["id"]: e for e in manifest["references"]}


def _seq(refs, name):
    header = next(h for h in refs if h.split("|")[0] == name)
    return refs[header].upper()


def test_all_four_references_are_present(resolver_refs_fasta):
    names = {h.split("|")[0] for h in resolver_refs_fasta}
    assert names == set(EXPECTED)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_reference_header_matches_expected_pair_and_position(resolver_refs_fasta, name):
    header = next(h for h in resolver_refs_fasta if h.split("|")[0] == name)
    meta = parse_resolver_meta(header)
    exp_pair, exp_pos, _, _ = EXPECTED[name]
    assert meta["pair"] == exp_pair
    assert meta["pos"] == exp_pos


# --- whole-sequence identity -------------------------------------------
#
# Regression: cps2K and cps14K shipped two bases short of their source
# records. Every position-level test still passed, because the deletion sat
# downstream of the diagnostic site. Only a whole-sequence check catches it.

@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_reference_sequence_matches_its_pinned_checksum(resolver_refs_fasta, manifest_by_id, name):
    seq = _seq(resolver_refs_fasta, name)
    digest = hashlib.sha256(seq.encode("ascii")).hexdigest()
    assert digest == manifest_by_id[name]["sha256"], (
        f"{name} is not the pinned sequence. Run "
        f"`python scripts/resolver_refs.py check` for detail.")


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_reference_length_matches_its_source_coordinates(resolver_refs_fasta, manifest_by_id, name):
    entry = manifest_by_id[name]
    seq = _seq(resolver_refs_fasta, name)
    span = entry["source_end"] - entry["source_start"] + 1
    assert len(seq) == entry["length"] == span, (
        f"{name}: {len(seq)} bases shipped, manifest says {entry['length']}, "
        f"{entry['accession']}:{entry['source_start']}-{entry['source_end']} spans {span}")


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_reference_is_an_intact_coding_sequence(resolver_refs_fasta, name):
    """A frameshifted reference misaligns the whole 3' end of the gene."""
    seq = _seq(resolver_refs_fasta, name)
    assert len(seq) % 3 == 0, f"{name}: length {len(seq)} is not a multiple of 3"
    assert seq[:3] == "ATG", f"{name}: starts {seq[:3]!r}, not ATG"
    assert seq[-3:] in STOPS, f"{name}: ends {seq[-3:]!r}, not a stop codon"
    internal = [i // 3 + 1 for i in range(0, len(seq) - 3, 3) if seq[i:i + 3] in STOPS]
    assert not internal, f"{name}: internal stop codon(s) at codon {internal}"


def test_shipped_fasta_passes_the_offline_manifest_check():
    """The same check the regeneration command runs, wired into the suite."""
    assert resolver_refs.check() == []


def test_manifest_records_provenance_for_every_reference(manifest):
    """accession.version, coordinates, strand, position, checksum, method."""
    assert manifest["extraction_method"]["command"]
    assert manifest["extraction_method"]["offline_check"]
    for entry in manifest["references"]:
        assert "." in entry["accession"], "accession must carry its version"
        assert entry["source_start"] < entry["source_end"]
        assert entry["strand"] in ("+", "-")
        assert entry["diagnostic_position"] > 0
        assert len(entry["sha256"]) == 64


def test_repair_is_a_no_op_on_the_shipped_references(tmp_path):
    """The two-base repairs are gated on the damaged checksum, so they can
    neither run twice nor touch an already-correct file.

    Run against a COPY: a repair() that found something to do would rewrite
    data/suis_resolver_refs.fasta and then report success, quietly healing the
    very state the rest of this file exists to detect.
    """
    import shutil
    original = Path(load_config()["resolver_refs_fasta"])
    copy = tmp_path / original.name
    shutil.copy(original, copy)

    assert resolver_refs.repair(copy) == []
    assert copy.read_bytes() == original.read_bytes()


def test_repair_restores_the_documented_damage_and_refuses_to_repeat(tmp_path):
    """The repair path itself, exercised on a deliberately damaged copy."""
    import shutil
    original = Path(load_config()["resolver_refs_fasta"])
    copy = tmp_path / original.name
    shutil.copy(original, copy)

    # Reintroduce the exact two-base deletions that shipped.
    records = resolver_refs.read_fasta(copy)
    damaged = []
    for header, seq in records:
        at = resolver_refs.REPAIRS.get(resolver_refs.ref_id(header), {}).get("insert_at")
        damaged.append((header, seq[:at] + seq[at + 2:]) if at else (header, seq))
    resolver_refs.write_fasta(copy, damaged)
    assert resolver_refs.check(copy) != [], "the damaged copy must fail the check"

    applied = resolver_refs.repair(copy)

    assert len(applied) == 2
    assert resolver_refs.check(copy) == []
    assert copy.read_bytes() == original.read_bytes(), "byte-for-byte restored"
    assert resolver_refs.repair(copy) == [], "and it will not run a second time"


# --- polarity and position ---------------------------------------------

@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_diagnostic_site_holds_the_expected_base(resolver_refs_fasta, name):
    seq = _seq(resolver_refs_fasta, name)
    _, pos, exp_codon, _ = EXPECTED[name]
    assert seq[pos - 1] == exp_codon[-1]


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_diagnostic_site_is_the_third_base_of_a_trp_or_cys_codon(resolver_refs_fasta, name):
    """The declared position must be the wobble base of codon 161.

    This is what makes the differing positions (483 vs 492) legitimate rather
    than a typo: the 1/14 references carry 9 extra bases at the 5' end.
    """
    seq = _seq(resolver_refs_fasta, name)
    _, pos, exp_codon, expected_serotype = EXPECTED[name]
    codon = seq[pos - 3:pos]
    assert codon == exp_codon
    assert codon in TRP | CYS, f"{name}: codon at {pos} is {codon!r}, not Trp/Cys"
    assert pos % 3 == 0, f"{name}: position {pos} is not a codon-third base"
    if codon in TRP:
        assert expected_serotype in ("2", "14"), "Trp161 -> Gal -> serotype 2 / 14"
    else:
        assert expected_serotype in ("1/2", "1"), "Cys161 -> GalNAc -> serotype 1/2 / 1"


def test_the_two_families_declare_different_positions_on_purpose(manifest_by_id):
    """483 vs 492 is a real 9 bp 5' offset, not a typo. Do not unify them."""
    assert manifest_by_id["cps2K"]["diagnostic_position"] == 483
    assert manifest_by_id["cps1/2K"]["diagnostic_position"] == 483
    assert manifest_by_id["cps14K"]["diagnostic_position"] == 492
    assert manifest_by_id["cps1L"]["diagnostic_position"] == 492
    assert manifest_by_id["cps14K"]["length"] - manifest_by_id["cps2K"]["length"] == 9


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_interpret_resolver_returns_the_published_serotype(resolver_refs_fasta, name):
    """Positive control, one per resolver output, read off the shipped file."""
    header = next(h for h in resolver_refs_fasta if h.split("|")[0] == name)
    seq = _seq(resolver_refs_fasta, name)
    _, pos, _, expected_serotype = EXPECTED[name]
    codon = seq[pos - 3:pos]
    ev = {"ref_id": header, "triplet": codon, "triplet_status": "OK",
          "base": codon[-1], "coding_status": "INTACT"}
    assert interpret_resolver(ev, {}) == expected_serotype


def test_the_two_pairs_share_flanking_context(resolver_refs_fasta):
    """Sanity-check the 9 bp offset: the sequence around the diagnostic base is
    the same in both pairs, which is why the position differs but the residue
    does not."""
    contexts = set()
    for name, (_, pos, _, _) in EXPECTED.items():
        seq = _seq(resolver_refs_fasta, name)
        contexts.add(seq[pos - 9:pos - 1] + seq[pos:pos + 8])
    assert len(contexts) == 1, f"flanking context diverges between references: {contexts}"
