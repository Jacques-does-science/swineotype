"""Real-BLAST regression tests over reference-derived fixtures.

These run blastn and makeblastdb for real. They are skipped -- loudly -- when
those tools are absent, because a mocked run of this code would not check the
thing under test: how the pipeline behaves on an actual alignment.

Scope: these fixtures are built from the project's own references, so they
establish correct behaviour on sequences of a given shape. They are not a
sample of any population and support no accuracy estimate.
"""
from __future__ import annotations

import pytest

from blastfixtures import (
    blast_config,
    build_assembly,
    cps2k_sequence,
    damaged_resolver_fasta,
    requires_blast,
    with_base_at,
)
from swineotype.main import process_one

pytestmark = requires_blast

# Sequence index 482, zero-based: the wobble base of the diagnostic codon.
DIAGNOSTIC_INDEX0 = 482


@pytest.fixture(scope="module")
def g_copy():
    """The corrected public cps2K sequence: TGG -> Trp161 -> serotype 2."""
    seq = cps2k_sequence()
    assert seq[480:483] == "TGG"
    return seq


@pytest.fixture(scope="module")
def t_copy(g_copy):
    """Identical except that sequence index 482 is T: TGT -> serotype 1/2."""
    seq = with_base_at(g_copy, DIAGNOSTIC_INDEX0, "T")
    assert seq[480:483] == "TGT"
    assert sum(a != b for a, b in zip(seq, g_copy)) == 1
    return seq


def run(tmp_path, name, copies, resolver_refs=None):
    asm = build_assembly(tmp_path / f"{name}.fasta", copies)
    cfg = blast_config(tmp_path / name)
    if resolver_refs is not None:
        cfg["resolver_refs_fasta"] = resolver_refs
    return process_one(str(asm), tmp_path / name / "out", 1, cfg)


# --- the conflict fixture ----------------------------------------------

@pytest.fixture(params=["corrected", "damaged"])
def reference_set(request, tmp_path_factory):
    """Both reference sets: the shipped one, and the one that shipped before.

    The conflict must be reported either way. At the reviewed baseline this
    same input returned a confident 1/2 against the damaged references and a
    confident 2 against the corrected ones, with no warning in either case --
    the answer depended on which reference happened to win on bit score.
    """
    if request.param == "corrected":
        return None
    return damaged_resolver_fasta(tmp_path_factory.mktemp("refs") / "damaged.fasta")


def test_two_disagreeing_copies_report_a_conflict(tmp_path, g_copy, t_copy, reference_set):
    row = run(tmp_path, "conflict", {"G_copy": g_copy, "T_copy": t_copy},
              resolver_refs=reference_set)

    assert row["stage2_status"] == "CONFLICTING_COPIES"
    assert row["final_serotype"] == "", "no exact label may be chosen between them"
    assert "conflicting_resolver_copies" in row["warnings"]
    # Both pieces of evidence survive, with their own coordinates.
    assert "G_copy" in row["resolver_loci"] and "T_copy" in row["resolver_loci"]
    assert "TGG" in row["resolver_loci"] and "TGT" in row["resolver_loci"]
    # The family-level result is still established and still reported.
    assert row["stage1_family"] == "2_vs_1_2"
    assert row["family_serotype"] == "2 or 1/2"
    assert row["status"] == "FAMILY_ONLY"


# --- single-copy controls ----------------------------------------------

def test_single_g_copy_calls_serotype_2(tmp_path, g_copy):
    row = run(tmp_path, "single_g", {"G_copy": g_copy})
    assert row["status"] == "STAGE2"
    assert row["final_serotype"] == "2"
    assert row["triplet"] == "TGG"
    assert row["triplet_status"] == "OK"
    assert row["warnings"] == ""


def test_single_t_copy_calls_serotype_1_2(tmp_path, t_copy):
    row = run(tmp_path, "single_t", {"T_copy": t_copy})
    assert row["status"] == "STAGE2"
    assert row["final_serotype"] == "1/2"
    assert row["triplet"] == "TGT"
    assert row["warnings"] == ""


def test_synonymous_copies_agree_on_serotype_1_2(tmp_path, g_copy):
    """Regression: a TGT copy and a TGC copy both encode Cys161 and both mean
    1/2, but the conflict check compared codon spellings and withheld the
    call. Built from the same corrected cps2K as the conflict fixture."""
    tgt = with_base_at(g_copy, DIAGNOSTIC_INDEX0, "T")
    tgc = with_base_at(g_copy, DIAGNOSTIC_INDEX0, "C")
    assert (tgt[480:483], tgc[480:483]) == ("TGT", "TGC")

    row = run(tmp_path, "synonymous", {"TGT_copy": tgt, "TGC_copy": tgc})

    assert row["stage2_status"] == "OK"
    assert row["final_serotype"] == "1/2"
    assert "conflicting" not in row["warnings"]
    assert "TGT_copy" in row["resolver_loci"] and "TGC_copy" in row["resolver_loci"], \
        "both copies are still reported"


def test_duplicate_equivalent_copies_do_not_conflict(tmp_path, g_copy):
    """Two identical copies agree, so there is nothing to withhold. Several
    references hitting them is not evidence of disagreement either."""
    row = run(tmp_path, "duplicate", {"G_copy": g_copy, "G_copy2": g_copy})
    assert row["stage2_status"] == "OK"
    assert row["final_serotype"] == "2"
    assert "conflicting" not in row["warnings"]


def test_several_references_hitting_one_locus_are_one_locus(tmp_path, g_copy):
    """cps2K and cps1/2K are 99% identical, so both align to the single copy.
    Different query names are not evidence of multiple physical copies."""
    row = run(tmp_path, "one_locus", {"G_copy": g_copy})
    assert row["resolver_loci"].count(";") == 0, row["resolver_loci"]
    assert row["final_serotype"] == "2"


# --- an undocumented codon ---------------------------------------------

def test_agg_at_the_diagnostic_codon_is_not_called(tmp_path, g_copy):
    """Regression: only the wobble base was read, so AGG passed as TGG.

    Sequence index 480 is the first base of the codon; making it A leaves the
    wobble base G untouched and turns Trp (TGG) into Arg (AGG).
    """
    agg = with_base_at(g_copy, 480, "A")
    assert agg[480:483] == "AGG"
    row = run(tmp_path, "agg", {"AGG_copy": agg})

    assert row["base"] == "G", "the wobble base alone still reads G"
    assert row["triplet"] == "AGG"
    assert row["triplet_status"] == "UNEXPECTED_CODON"
    assert row["final_serotype"] == "", "an Arg codon is not a documented state"
    assert row["stage2_status"] == "INVALID_TRIPLET"
    assert row["family_serotype"] == "2 or 1/2", "the family result survives"


def test_deleted_diagnostic_codon_is_not_called(tmp_path, g_copy):
    deleted = g_copy[:480] + g_copy[483:]
    row = run(tmp_path, "deleted", {"DEL_copy": deleted})
    assert row["final_serotype"] == ""
    assert row["stage2_status"] in ("INVALID_TRIPLET", "CODING_DISRUPTED")


def test_no_cps_match_reports_no_species(tmp_path):
    """An assembly with no cps content at all. Regression: this returned
    `species=Streptococcus suis` because the column defaulted to the target."""
    from blastfixtures import write_fasta
    asm = write_fasta(tmp_path / "empty.fasta", {"c1": "ACGT" * 500})
    cfg = blast_config(tmp_path / "nohit")
    row = process_one(str(asm), tmp_path / "nohit" / "out", 1, cfg)

    assert row["status"] == "NO_CPS_MATCH"
    assert row["species"] == ""
    assert row["matched_reference_taxon"] == ""
    assert row["final_serotype"] == ""


# --- positive controls, one per resolver output ------------------------

RESOLVER_OUTPUTS = [
    # (resolver reference, wzx/wzy source record, expected serotype)
    ("cps2K", "BR001000", "2"),
    ("cps1/2K", "AB737816", "1/2"),
    ("cps14K", "AB737822", "14"),
    ("cps1L", "AB737817", "1"),
]


@pytest.mark.parametrize("ref_name, marker_accession, expected", RESOLVER_OUTPUTS)
def test_each_resolver_output_is_reproduced_end_to_end(tmp_path, ref_name,
                                                       marker_accession, expected):
    """One positive control per documented output of the resolver.

    The assembly is the shipped reference itself: its own cpsK/cps1L gene plus
    the wzx/wzy markers of the same source record. Anything other than
    `expected` means the reference, the family selection or the codon
    interpretation has moved.
    """
    from blastfixtures import by_id, markers_from, resolver_records

    seq = by_id(resolver_records())[ref_name]
    asm = build_assembly(tmp_path / "control.fasta", {"cpsK_copy": seq},
                         markers=markers_from(marker_accession))
    row = process_one(str(asm), tmp_path / "out", 1, blast_config(tmp_path))

    assert row["status"] == "STAGE2", row
    assert row["final_serotype"] == expected
    assert row["triplet_status"] == "OK"
    assert row["triplet"] in ("TGG", "TGT", "TGC")
