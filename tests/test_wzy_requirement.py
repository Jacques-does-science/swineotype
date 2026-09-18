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

from swineotype.config import load_config
from swineotype.stages import parse_whitelist_headers, stage1_score

SUIS = "Streptococcus suis"


def cfg(**over):
    c = {"min_pid": 85.0, "min_cov": 0.8, "plurality": 0.6, "delta": 100,
         "ambig_set": {"1", "14", "2", "1/2"}, "keep_debug": False,
         "tmp_dir": "tmp", "require_wzy": 1}
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

def make_main_config(**over):
    c = {"target_species": SUIS, "pair_1_14": {"1", "14"}, "pair_2_1_2": {"2", "1/2"},
         "ambig_set": {"1", "14", "2", "1/2"}, "delta": 100, "require_wzy": 1,
         "wzxwzy_fasta": "w.fasta", "resolver_refs_fasta": "r.fasta"}
    c.update(over)
    return c


def test_no_wzy_match_is_reported_as_a_candidate_novel_locus(monkeypatch, tmp_path):
    from swineotype import main as main_mod
    s1 = {"top": None, "second": None, "decisive": False, "delta": 0.0,
          "must_stage2_for_pair": False, "top_species": None, "top_wzx_only": "27"}
    monkeypatch.setattr(main_mod, "stage1_score", lambda *a, **k: s1)
    monkeypatch.setattr(main_mod, "ensure_unix_line_endings", lambda p, t: p)
    monkeypatch.setattr(main_mod, "stage2_resolver_call",
                        lambda *a, **k: pytest.fail("Stage 2 must not run without a callable type"))

    r = main_mod.process_one("iso.fasta", tmp_path, 1, make_main_config(tmp_dir=tmp_path))

    assert r["status"] == "NO_WZY_MATCH"
    assert r["final_serotype"] == ""
    assert r["stage1_top"] == "", "no serotype may be implied"
    assert r["stage1_wzx_only"] == "27", "the lead is reported separately"
    assert r["species"] == "", "wzx is not species-discriminating either"
    assert "candidate_novel_capsular_locus" in r["warnings"]
    assert "nearest_wzx_relative:cps_type_27" in r["warnings"]


def test_stage2_status_is_surfaced(monkeypatch, tmp_path):
    """It was computed in process_one and then discarded, so a user could not
    see why Stage 2 produced nothing."""
    from swineotype import main as main_mod
    s1 = {"top": "2", "second": "1/2", "decisive": False, "delta": 10.0,
          "must_stage2_for_pair": True, "top_species": SUIS, "top_wzx_only": None}
    monkeypatch.setattr(main_mod, "stage1_score", lambda *a, **k: s1)
    monkeypatch.setattr(main_mod, "ensure_unix_line_endings", lambda p, t: p)
    monkeypatch.setattr(main_mod, "stage2_resolver_call", lambda *a, **k: None)

    r = main_mod.process_one("iso.fasta", tmp_path, 1, make_main_config(tmp_dir=tmp_path))
    assert r["stage2_status"] == "NO_HSP_OR_LOW_QUAL"
    assert r["status"] == "NO_CALL_STAGE2"


def test_stage2_status_ok_on_a_successful_resolve(monkeypatch, tmp_path):
    from swineotype import main as main_mod
    s1 = {"top": "2", "second": "1/2", "decisive": False, "delta": 10.0,
          "must_stage2_for_pair": True, "top_species": SUIS, "top_wzx_only": None}
    ev = {"ref_id": "cps2K|pair=2_vs_1_2|pos=483|G_serotype=2|CT_serotype=1/2",
          "contig": "c1", "contig_pos": 883, "strand": "+", "base": "G"}
    monkeypatch.setattr(main_mod, "stage1_score", lambda *a, **k: s1)
    monkeypatch.setattr(main_mod, "ensure_unix_line_endings", lambda p, t: p)
    monkeypatch.setattr(main_mod, "stage2_resolver_call", lambda *a, **k: ev)

    r = main_mod.process_one("iso.fasta", tmp_path, 1, make_main_config(tmp_dir=tmp_path))
    assert r["stage2_status"] == "OK"
    assert r["final_serotype"] == "2"
