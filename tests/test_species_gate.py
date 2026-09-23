"""Species gate: a cps type belonging to another organism is not a serotype.

Six serotypes of the original 35 have been moved out of S. suis. Their cps
loci are still needed in the panel -- they are what a hit to that organism
looks like -- but reporting one as "S. suis serotype 22" is wrong at the
species level. These tests pin both the reference labels and the gate.

  20, 22, 26 -> Streptococcus parasuis     (Nomoto et al. 2015, IJSEM 65:438)
  33         -> Streptococcus ruminantium  (Tohya et al. 2017, IJSEM 67:3660)
  32, 34     -> Streptococcus orisratti    (Hill et al. 2005, Vet Microbiol 107:63)

Leaving 29 true S. suis serotypes: 1-19, 21, 23-25, 27-31 and 1/2.
"""
import pytest

from helpers import SUIS, main_config, stage1_double
from swineotype.config import load_config
from swineotype.stages import header_tag, parse_whitelist_headers

REASSIGNED = {
    "20": "Streptococcus parasuis",
    "22": "Streptococcus parasuis",
    "26": "Streptococcus parasuis",
    "32": "Streptococcus orisratti",
    "33": "Streptococcus ruminantium",
    "34": "Streptococcus orisratti",
}

TRUE_SUIS_SEROTYPES = {str(n) for n in list(range(1, 20)) + [21, 23, 24, 25, 27, 28, 29, 30, 31]} | {"1/2"}


@pytest.fixture(scope="module")
def type_to_species():
    _, _, mapping = parse_whitelist_headers(load_config()["wzxwzy_fasta"])
    return mapping


# --- header parsing -----------------------------------------------------

def test_header_tag_reads_values_containing_spaces():
    """Species binomials and literature refs contain spaces, so the tag
    parser cannot be a whitespace split."""
    h = "wzy_X [type_id=22] [species=Streptococcus parasuis] [reassignment_ref=Nomoto et al. 2015, IJSEM 65:438]"
    assert header_tag(h, "species") == "Streptococcus parasuis"
    assert header_tag(h, "type_id") == "22"
    assert header_tag(h, "reassignment_ref") == "Nomoto et al. 2015, IJSEM 65:438"
    assert header_tag(h, "absent") is None


def test_header_tag_does_not_match_a_prefix():
    h = "wzy_X [type_id=22] [type_id_extra=99]"
    assert header_tag(h, "type_id") == "22"


# --- the shipped reference panel ---------------------------------------

def test_every_reference_declares_a_species(type_to_species):
    missing = [t for t, sp in type_to_species.items() if not sp]
    assert not missing, f"cps types with no [species=] tag: {missing}"


@pytest.mark.parametrize("cps_type, species", sorted(REASSIGNED.items()))
def test_reassigned_types_are_labelled_with_their_current_species(type_to_species, cps_type, species):
    assert type_to_species[cps_type] == species


def test_exactly_the_29_true_serotypes_are_labelled_s_suis(type_to_species):
    labelled_suis = {t for t, sp in type_to_species.items() if sp == SUIS}
    assert labelled_suis == TRUE_SUIS_SEROTYPES
    assert len(labelled_suis) == 29


def test_reassigned_types_are_kept_not_deleted(type_to_species):
    """They must stay in the panel: they are how the tool recognises those
    organisms. Deleting them would turn a wrong answer into a silent NO_CALL."""
    assert set(REASSIGNED) <= set(type_to_species)


# --- the gate -----------------------------------------------------------

@pytest.mark.parametrize("cps_type, species", sorted(REASSIGNED.items()))
def test_gate_rejects_a_non_target_species(patched_stages, tmp_path, cps_type, species):
    main_mod = patched_stages(stage1_double(top=cps_type, top_species=species),
                              forbid_stage2=True)

    row = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert row["status"] == "NON_TARGET_SPECIES"
    assert row["matched_reference_taxon"] == species
    assert row["species"] == species, "deprecated alias still carries it"
    assert row["final_serotype"] == "", "must not report a serotype for another organism"
    assert row["stage1_top"] == cps_type, "the cps type is still reported, as evidence"
    assert species.replace(" ", "_") in row["warnings"]


def test_gate_lets_s_suis_through(patched_stages, tmp_path):
    main_mod = patched_stages(stage1_double(top="9", second="7"))

    row = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert row["status"] == "STAGE1"
    assert row["final_serotype"] == "9"
    assert row["matched_reference_taxon"] == SUIS


def test_target_species_is_configurable(patched_stages, tmp_path):
    """The gate compares against config, not a hard-coded string, so the panel
    can be repurposed without editing Python."""
    main_mod = patched_stages(stage1_double(top="22", top_species="Streptococcus parasuis"))

    cfg = main_config(target_species="Streptococcus parasuis", tmp_dir=tmp_path)
    row = main_mod.process_one("iso.fasta", tmp_path, 1, cfg)

    assert row["status"] == "STAGE1"
    assert row["final_serotype"] == "22"
    assert row["matched_reference_taxon"] == "Streptococcus parasuis"


def test_unlabelled_type_is_not_blocked(patched_stages, tmp_path):
    """A reference with no [species=] tag should not be treated as foreign --
    fail open on the serotype, since the panel is S. suis by construction."""
    main_mod = patched_stages(stage1_double(top="9", top_species=None))

    row = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))
    assert row["status"] == "STAGE1"
    assert row["final_serotype"] == "9"
    assert row["matched_reference_taxon"] == ""


# --- the matched taxon is not an identification of the input ------------

def test_matched_reference_taxon_is_not_reported_as_the_input_species(patched_stages, tmp_path):
    """Reference metadata says what the reference is, not what the input is."""
    main_mod = patched_stages(stage1_double(top="9"))

    row = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path))

    assert row["matched_reference_taxon"] == SUIS
    assert row["input_species"] == "", "nothing measured the input"
    assert row["species_assessment"] == "NOT_ASSESSED"


def test_supplied_input_species_is_recorded_separately(patched_stages, tmp_path):
    main_mod = patched_stages(stage1_double(top="9"))

    row = main_mod.process_one("iso.fasta", tmp_path, 1, main_config(tmp_dir=tmp_path),
                               input_species="Streptococcus suis")

    assert row["input_species"] == "Streptococcus suis"
    assert row["species_assessment"] == "USER_SUPPLIED"
    assert row["matched_reference_taxon"] == SUIS, "still a separate column"
