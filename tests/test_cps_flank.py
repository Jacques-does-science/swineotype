"""Telling "this cps region cannot be typed" from "there is no cps region".

Both used to come back as NO_CPS_MATCH, so a Streptococcus suis genome whose
cassette is absent from the panel -- a novel cps locus -- was indistinguishable
from an unrelated organism or an unusable assembly. The conserved flanking
region separates them.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from helpers import main_config, stage1_double

NO_MATCH = dict(top=None, family_top=None, family_decisive=False,
                top_species=None, top_wzx_only=None)


def flank(present, *, found=(), span=None):
    return {"present": present, "flanks_found": list(found),
            "coverage": {}, "interflank_span": span}


# --- which status, and on what evidence --------------------------------

def test_flanks_present_is_not_the_same_as_nothing_found(patched_stages, tmp_path):
    main_mod = patched_stages(stage1_double(**NO_MATCH), forbid_stage2=True,
                              flank=flank(True, found=["cps_flank_5prime",
                                                       "cps_flank_3prime"],
                                          span=11703))
    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "CPS_REGION_UNTYPED"
    assert "no_cps_reference_matched" in r["warnings"]
    assert "cps_flanking_region_present:cps_flank_5prime+cps_flank_3prime" in r["warnings"]
    assert "interflank_span=11703" in r["warnings"]
    assert "cps_cassette_not_matched_by_panel" in r["warnings"]


def test_nothing_at_all_still_reports_no_cps_match(patched_stages, tmp_path):
    """The existing status keeps its meaning, narrowed: no type AND no flank."""
    main_mod = patched_stages(stage1_double(**NO_MATCH), forbid_stage2=True,
                              flank=flank(False))
    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "NO_CPS_MATCH"
    assert "no_cps_flanking_region" in r["warnings"]


def test_an_unmeasurable_span_is_said_to_be_unmeasured(patched_stages, tmp_path):
    """On a draft assembly the two flanks usually land on different contigs.
    Silence there would read as "no cassette", which is not what was found."""
    main_mod = patched_stages(stage1_double(**NO_MATCH), forbid_stage2=True,
                              flank=flank(True, found=["cps_flank_5prime"]))
    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "CPS_REGION_UNTYPED"
    assert "interflank_span_not_assessed" in r["warnings"]
    assert "interflank_span=" not in r["warnings"]


def test_neither_status_claims_a_species_or_a_serotype(patched_stages, tmp_path):
    """A flank is chromosomal context, not an identification, and certainly
    not a type. Nothing here may populate the species or serotype columns."""
    for f in (flank(True, found=["cps_flank_5prime", "cps_flank_3prime"]),
              flank(False)):
        main_mod = patched_stages(stage1_double(**NO_MATCH), forbid_stage2=True,
                                  flank=f)
        r = main_mod.process_one("iso.fasta", tmp_path, 1,
                                 main_config(tmp_dir=tmp_path))
        assert r["final_serotype"] == ""
        assert r["matched_reference_taxon"] == ""
        assert r["species"] == ""
        assert r["species_assessment"] == "NOT_ASSESSED"


def test_novelty_is_never_asserted(patched_stages, tmp_path):
    """The flanks survive a deleted cassette, so their presence cannot
    establish that a novel one is there. Report evidence, not a conclusion."""
    main_mod = patched_stages(stage1_double(**NO_MATCH), forbid_stage2=True,
                              flank=flank(True, found=["cps_flank_5prime",
                                                       "cps_flank_3prime"],
                                          span=11703))
    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert "novel" not in r["warnings"].lower()


def test_the_probe_does_not_run_when_stage1_matched(patched_stages, monkeypatch,
                                                    tmp_path):
    """It is a fallback, not a per-sample cost on the common path."""
    main_mod = patched_stages(stage1_double(top="9", family_top="type:9"),
                              forbid_stage2=True)
    monkeypatch.setattr(main_mod, "cps_flank_probe", lambda *a, **k: pytest.fail(
        "the flank probe must not run once a type matched"))

    r = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert r["status"] == "STAGE1"
    assert r["final_serotype"] == "9"


# --- the probe itself, against real sequence ---------------------------

from blastfixtures import blast_config, requires_blast, write_fasta  # noqa: E402

from swineotype.stages import cps_flank_probe  # noqa: E402

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


@requires_blast
def test_the_probe_finds_both_flanks_in_a_real_suis_genome(tmp_path):
    genome = EXAMPLES / "Streptococcus" / "GCA_000014325.1_ASM1432v1_genomic.fasta"
    if not genome.exists():
        pytest.skip(f"example genome not shipped: {genome}")
    cfg = blast_config(tmp_path)
    out = tmp_path / "run"; out.mkdir()
    r = cps_flank_probe(str(genome), str(cfg["cps_flank_fasta"]), 1, out, cfg)

    assert r["present"] is True
    assert set(r["flanks_found"]) == {"cps_flank_5prime", "cps_flank_3prime"}
    # A complete genome is one contig, so the cassette between them is
    # measurable and is of the order a cps locus actually is.
    assert r["interflank_span"] is not None
    assert 10_000 < r["interflank_span"] < 60_000


@requires_blast
def test_the_probe_finds_nothing_in_an_unrelated_organism(tmp_path):
    genome = EXAMPLES / "Actinobacillus" / "CP031873.1.fasta"
    if not genome.exists():
        pytest.skip(f"example genome not shipped: {genome}")
    cfg = blast_config(tmp_path)
    out = tmp_path / "run"; out.mkdir()
    r = cps_flank_probe(str(genome), str(cfg["cps_flank_fasta"]), 1, out, cfg)

    assert r["present"] is False
    assert r["flanks_found"] == []


@requires_blast
def test_the_probe_finds_nothing_in_junk(tmp_path):
    asm = write_fasta(tmp_path / "junk.fasta", {"c1": "ACGT" * 500})
    cfg = blast_config(tmp_path)
    out = tmp_path / "run"; out.mkdir()
    r = cps_flank_probe(str(asm), str(cfg["cps_flank_fasta"]), 1, out, cfg)

    assert r["present"] is False
