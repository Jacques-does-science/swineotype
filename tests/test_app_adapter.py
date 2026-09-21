"""APP adapter: identities, staged filenames, and reading the suis summary."""
import csv
from pathlib import Path

import pandas as pd
import pytest

from swineotype.adapters.app import (
    STAGED_SUFFIX,
    read_swineotype_summary,
    stage_assembly,
    unique_sample_names,
)
from swineotype.main import SUMMARY_COLUMNS


# --- reading the swineotype summary ------------------------------------

def write_summary(path: Path, rows, delimiter=","):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS,
                                delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({**{c: "" for c in SUMMARY_COLUMNS}, **r})


def test_a_csv_summary_is_read_as_csv(tmp_path):
    """Regression: `pd.read_csv(path, sep="\\t")` on a CSV does not raise -- it
    returns one column whose name is the whole header line. The try/except
    fallback therefore never ran, and the merge died with KeyError('sample')."""
    p = tmp_path / "summary.csv"
    write_summary(p, [{"sample": "iso1", "final_serotype": "2"}])

    df = read_swineotype_summary(p)

    assert "sample" in df.columns, "the merge key must survive the read"
    assert df.loc[0, "sample"] == "iso1"
    assert len(df.columns) == len(SUMMARY_COLUMNS)


def test_merging_on_sample_works_after_the_fix(tmp_path):
    """The failure this actually caused, end to end."""
    p = tmp_path / "summary.csv"
    write_summary(p, [{"sample": "iso1"}, {"sample": "iso2"}])
    app = pd.DataFrame({"sample": ["iso1", "iso3"], "app_serovar": ["APP_5", "APP_7"]})

    merged = read_swineotype_summary(p).merge(app, on="sample", how="outer")

    assert set(merged["sample"]) == {"iso1", "iso2", "iso3"}
    assert merged.set_index("sample").loc["iso1", "app_serovar"] == "APP_5"


def test_a_tab_separated_summary_is_still_read(tmp_path):
    p = tmp_path / "summary.tsv"
    write_summary(p, [{"sample": "iso1"}], delimiter="\t")
    assert read_swineotype_summary(p).loc[0, "sample"] == "iso1"


def test_a_file_without_a_sample_column_is_rejected_by_name(tmp_path):
    p = tmp_path / "other.csv"
    p.write_text("name,value\niso1,2\n")
    with pytest.raises(ValueError, match="no 'sample' column"):
        read_swineotype_summary(p)


def test_an_empty_summary_is_rejected(tmp_path):
    p = tmp_path / "empty.csv"
    p.write_text("")
    with pytest.raises(ValueError, match="empty"):
        read_swineotype_summary(p)


def test_numeric_looking_sample_names_stay_strings(tmp_path):
    """A merge key that pandas reads as int64 on one side and object on the
    other silently matches nothing."""
    p = tmp_path / "summary.csv"
    write_summary(p, [{"sample": "0012"}])
    assert read_swineotype_summary(p).loc[0, "sample"] == "0012"


# --- sample identity ---------------------------------------------------

def test_two_assembly_fasta_paths_get_distinct_identities(tmp_path):
    """Regression: both reduced to the stem 'assembly', so they shared one row
    in the sample sheet and one staged symlink."""
    a = tmp_path / "runA" / "assembly.fasta"
    b = tmp_path / "runB" / "assembly.fasta"
    for p in (a, b):
        p.parent.mkdir()
        p.write_text(">c\nACGT\n")

    names = unique_sample_names([a, b])

    assert len(set(names)) == 2
    assert all(n.startswith("assembly__") for n in names)


def test_unique_basenames_are_left_alone(tmp_path):
    a = tmp_path / "iso1.fasta"; a.write_text(">c\nACGT\n")
    b = tmp_path / "iso2.fna"; b.write_text(">c\nACGT\n")
    assert unique_sample_names([a, b]) == ["iso1", "iso2"]


def test_identity_is_stable_for_the_same_path(tmp_path):
    a = tmp_path / "x" / "assembly.fasta"; a.parent.mkdir(); a.write_text(">c\nA\n")
    b = tmp_path / "y" / "assembly.fasta"; b.parent.mkdir(); b.write_text(">c\nA\n")
    assert unique_sample_names([a, b]) == unique_sample_names([a, b])


# --- staging -----------------------------------------------------------

def test_a_stale_valid_symlink_is_replaced(tmp_path):
    """Regression: rerunning with a new input of the same name kept the old
    link, because it still resolved and `exists()` therefore said True."""
    old = tmp_path / "old.fasta"; old.write_text(">old\nA\n")
    new = tmp_path / "new.fasta"; new.write_text(">new\nC\n")
    dest = tmp_path / "staged.fasta"

    stage_assembly(old, dest)
    stage_assembly(new, dest)

    assert dest.resolve() == new.resolve()
    assert dest.read_text() == ">new\nC\n"


def test_a_broken_symlink_is_replaced_without_raising(tmp_path):
    """Regression: `exists()` follows symlinks, so a dangling link answered
    False and the symlink_to() that followed raised FileExistsError."""
    gone = tmp_path / "gone.fasta"; gone.write_text(">g\nA\n")
    dest = tmp_path / "staged.fasta"
    stage_assembly(gone, dest)
    gone.unlink()
    assert dest.is_symlink() and not dest.exists()

    real = tmp_path / "real.fasta"; real.write_text(">r\nC\n")
    stage_assembly(real, dest)

    assert dest.read_text() == ">r\nC\n"


def test_staging_is_idempotent(tmp_path):
    src = tmp_path / "a.fasta"; src.write_text(">a\nA\n")
    dest = tmp_path / "staged.fasta"
    stage_assembly(src, dest)
    stage_assembly(src, dest)
    assert dest.resolve() == src.resolve()


def test_staged_name_always_carries_the_fasta_suffix():
    """The workflow's input pattern is the literal `{sample}.fasta`, so a .fa
    or .fna staged under its own name left the rule with no input at all."""
    assert STAGED_SUFFIX == ".fasta"


def test_staged_filenames_are_normalised(tmp_path, monkeypatch):
    """End to end through the adapter's staging step, for .fa / .fna inputs."""
    from swineotype.adapters import app as app_mod

    inputs = [tmp_path / "iso1.fa", tmp_path / "iso2.fna", tmp_path / "iso3.fasta"]
    for p in inputs:
        p.write_text(">c\nACGT\n")
    staged_dir = tmp_path / "staged"; staged_dir.mkdir()

    names = app_mod.unique_sample_names(inputs)
    for src, name in zip(inputs, names):
        app_mod.stage_assembly(src, staged_dir / f"{name}{app_mod.STAGED_SUFFIX}")

    assert sorted(p.name for p in staged_dir.iterdir()) == \
        ["iso1.fasta", "iso2.fasta", "iso3.fasta"]


def test_sample_sheet_names_match_the_staged_filenames(tmp_path):
    """The sheet's sample_name is the wildcard the workflow substitutes into
    `{tmpdir}/assemblies/{sample}.fasta`; if they diverge, nothing runs."""
    inputs = [tmp_path / "runA" / "assembly.fa", tmp_path / "runB" / "assembly.fna"]
    for p in inputs:
        p.parent.mkdir(); p.write_text(">c\nACGT\n")

    names = unique_sample_names(inputs)
    staged = [f"{n}{STAGED_SUFFIX}" for n in names]

    assert [Path(s).stem for s in staged] == names
    assert len(set(staged)) == 2


# --- the adapter's setup phase, end to end -----------------------------

def fake_kma_db(third_party: Path):
    db = third_party / "db"
    db.mkdir(parents=True, exist_ok=True)
    for suffix in (".fasta", ".seq.b", ".comp.b", ".length.b"):
        (db / f"Actinobacillus_pleuropneumoniae{suffix}").write_text("x")
    cfg = third_party / "config"
    cfg.mkdir(parents=True, exist_ok=True)
    (cfg / "serovar_profiles.yaml").write_text("profiles: {}\n")


@pytest.fixture
def app_setup(tmp_path, monkeypatch):
    """Run run_app_analysis up to (not including) the Snakemake invocation."""
    from swineotype.adapters import app as app_mod

    third_party = tmp_path / "third_party" / "serovar_detector"
    fake_kma_db(third_party)
    monkeypatch.setattr(app_mod, "__file__",
                        str(tmp_path / "swineotype" / "adapters" / "app.py"))

    calls = {}

    class Done:
        returncode = 1  # stop before the results step

    def fake_run(cmd, *a, **k):
        calls["cmd"] = cmd
        return Done()

    monkeypatch.setattr(app_mod.subprocess, "run", fake_run)

    def go(patterns, out_dir):
        with pytest.raises(SystemExit):
            app_mod.run_app_analysis(assembly=patterns, out_dir=str(out_dir),
                                     threads=1, swineotype_summary=None)
        return calls

    return go


def test_setup_stages_every_input_as_sample_dot_fasta(tmp_path, app_setup):
    inputs = tmp_path / "in"; inputs.mkdir()
    for name in ("iso1.fa", "iso2.fna", "iso3.fasta"):
        (inputs / name).write_text(">c\nACGT\n")

    out = tmp_path / "out"
    app_setup([str(inputs / "*")], out)

    staged = out / "app_detector" / "tmp" / "assemblies"
    assert sorted(p.name for p in staged.iterdir()) == \
        ["iso1.fasta", "iso2.fasta", "iso3.fasta"]

    sheet = list(csv.DictReader((out / "app_detector" / "results" / "schemas"
                                 / "sample_sheet.csv").open()))
    assert [r["sample_name"] for r in sheet] == ["iso1", "iso2", "iso3"]
    assert {f"{r['sample_name']}.fasta" for r in sheet} == \
        {p.name for p in staged.iterdir()}, "sheet names must match staged filenames"


def test_setup_disambiguates_two_assembly_fasta_inputs(tmp_path, app_setup):
    for run_name in ("runA", "runB"):
        d = tmp_path / "in" / run_name; d.mkdir(parents=True)
        (d / "assembly.fasta").write_text(f">{run_name}\nACGT\n")

    out = tmp_path / "out"
    app_setup([str(tmp_path / "in" / "*" / "assembly.fasta")], out)

    staged = sorted(p.name for p in
                    (out / "app_detector" / "tmp" / "assemblies").iterdir())
    assert len(staged) == 2, "two different files must not share one staged name"


def test_setup_deduplicates_overlapping_patterns(tmp_path, app_setup):
    inputs = tmp_path / "in"; inputs.mkdir()
    (inputs / "iso1.fasta").write_text(">c\nACGT\n")

    out = tmp_path / "out"
    app_setup([str(inputs / "*.fasta"), str(inputs / "iso1.fasta")], out)

    sheet = list(csv.DictReader((out / "app_detector" / "results" / "schemas"
                                 / "sample_sheet.csv").open()))
    assert [r["sample_name"] for r in sheet] == ["iso1"], "one file, one row"


def test_setup_rejects_a_pattern_that_matched_nothing(tmp_path, monkeypatch, capsys):
    from swineotype.adapters import app as app_mod
    inputs = tmp_path / "in"; inputs.mkdir()
    (inputs / "iso1.fasta").write_text(">c\nACGT\n")

    with pytest.raises(SystemExit) as exc:
        app_mod.run_app_analysis(
            assembly=[str(inputs / "*.fasta"), str(inputs / "*.fna")],
            out_dir=str(tmp_path / "out"), threads=1, swineotype_summary=None)
    assert exc.value.code == 2


def test_setup_writes_a_run_record(tmp_path, app_setup):
    import json
    inputs = tmp_path / "in"; inputs.mkdir()
    (inputs / "iso1.fasta").write_text(">c\nACGT\n")

    out = tmp_path / "out"
    app_setup([str(inputs / "*.fasta")], out)

    record = json.loads((out / "swineotype_app_run.json").read_text())
    assert record["swineotype_version"]
    assert record["adapter"] == "serovar_detector"
    assert list(record["samples"]) == ["iso1"]
    assert record["workflow_config"]["threshold"] == 98.0
