"""APP adapter: identities, staging for serovar_detector, and merging with the suis summary."""
import csv
import json
import subprocess
from pathlib import Path

import pytest

from swineotype.adapters.app import merge_app_results, read_swineotype_summary
from swineotype.main import SUMMARY_COLUMNS
from swineotype.utils import unique_run_names


# --- the swineotype summary, and merging APP calls into it --------------

def write_summary(path: Path, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({**{c: "" for c in SUMMARY_COLUMNS}, **r})


def write_serovars_tsv(path: Path, calls):
    """serovars.tsv as serovar_detector 1.1.x writes it for assemblies."""
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Sample", "Mapper", "Suggested_serovar", "Frequency",
                         "Serovar_match", "Serovar_partial"])
        for sample, serovar in calls:
            writer.writerow([sample, "Blastn", serovar, "1 of 1", "cps", ""])
    return path


def read_rows(path: Path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def test_a_csv_summary_is_read_as_csv(tmp_path):
    """Regression: `pd.read_csv(path, sep="\\t")` on a CSV does not raise -- it
    returns one column whose name is the whole header line. The try/except
    fallback therefore never ran, and the merge died with KeyError('sample')."""
    p = tmp_path / "summary.csv"
    write_summary(p, [{"sample": "iso1", "final_serotype": "2"}])

    columns, rows = read_swineotype_summary(p)

    assert "sample" in columns, "the merge key must survive the read"
    assert columns == SUMMARY_COLUMNS
    assert rows[0]["sample"] == "iso1"


def test_app_calls_are_outer_joined_onto_the_summary(tmp_path):
    """The failure this actually caused, end to end through the adapter."""
    summary = tmp_path / "summary.csv"
    write_summary(summary, [{"sample": "iso1", "final_serotype": "2"},
                            {"sample": "iso2", "final_serotype": "1/2"}])
    tsv = write_serovars_tsv(tmp_path / "serovars.tsv", [("iso1", "S5"), ("iso3", "S7")])

    merge_app_results(summary, tsv)

    rows = read_rows(summary)
    assert [r["sample"] for r in rows] == ["iso1", "iso2", "iso3"], "ordered by sample"
    assert {r["sample"]: r["app_serovar"] for r in rows} == \
        {"iso1": "S5", "iso2": "", "iso3": "S7"}
    assert list(rows[0]) == SUMMARY_COLUMNS + ["app_serovar"]


def test_merging_keeps_serotypes_as_text(tmp_path):
    """Regression: the pandas merge re-typed all-numeric columns, so an outer
    join that introduced a blank row turned final_serotype "2" into "2.0"."""
    summary = tmp_path / "summary.csv"
    write_summary(summary, [{"sample": "iso1", "final_serotype": "2"},
                            {"sample": "iso2", "final_serotype": "14"}])
    tsv = write_serovars_tsv(tmp_path / "serovars.tsv", [("app_only", "S7")])

    merge_app_results(summary, tsv)

    finals = {r["sample"]: r["final_serotype"] for r in read_rows(summary)}
    assert finals == {"app_only": "", "iso1": "2", "iso2": "14"}


def test_without_a_summary_the_app_calls_are_written_in_their_own_order(tmp_path):
    summary = tmp_path / "summary.csv"
    tsv = write_serovars_tsv(tmp_path / "serovars.tsv", [("iso9", "S1"), ("iso2", "S3")])

    merge_app_results(summary, tsv)

    assert [(r["sample"], r["app_serovar"]) for r in read_rows(summary)] == \
        [("iso9", "S1"), ("iso2", "S3")]


def test_a_serovars_tsv_missing_its_columns_is_rejected(tmp_path):
    tsv = tmp_path / "serovars.tsv"
    tsv.write_text("Sample\tMapper\niso1\tBlastn\n")
    with pytest.raises(ValueError, match="Suggested_serovar"):
        merge_app_results(tmp_path / "summary.csv", tsv)


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
    """A merge key read as int64 on one side and text on the other matches
    nothing; "0012" and "12" are different samples."""
    summary = tmp_path / "summary.csv"
    write_summary(summary, [{"sample": "0012"}, {"sample": "12"}])
    tsv = write_serovars_tsv(tmp_path / "serovars.tsv", [("0012", "S2")])

    merge_app_results(summary, tsv)

    assert {r["sample"]: r["app_serovar"] for r in read_rows(summary)} == \
        {"0012": "S2", "12": ""}


# --- sample identity ---------------------------------------------------

def test_two_assembly_fasta_paths_get_distinct_identities(tmp_path):
    """Regression: both reduced to the stem 'assembly', so they shared one
    staged file and one sample identity."""
    a = tmp_path / "runA" / "assembly.fasta"
    b = tmp_path / "runB" / "assembly.fasta"
    for p in (a, b):
        p.parent.mkdir()
        p.write_text(">c\nACGT\n")

    names = unique_run_names([a, b])

    assert len(set(names)) == 2
    assert all(n.startswith("assembly__") for n in names)


def test_unique_basenames_are_left_alone(tmp_path):
    a = tmp_path / "iso1.fasta"; a.write_text(">c\nACGT\n")
    b = tmp_path / "iso2.fna"; b.write_text(">c\nACGT\n")
    assert unique_run_names([a, b]) == ["iso1", "iso2"]


def test_identity_is_stable_for_the_same_path(tmp_path):
    a = tmp_path / "x" / "assembly.fasta"; a.parent.mkdir(); a.write_text(">c\nA\n")
    b = tmp_path / "y" / "assembly.fasta"; b.parent.mkdir(); b.write_text(">c\nA\n")
    assert unique_run_names([a, b]) == unique_run_names([a, b])


# --- running serovar_detector ------------------------------------------

@pytest.fixture
def serovar_detector(monkeypatch):
    """Stand in for the serovar_detector command.

    Records each call and, like the real tool, writes <-o>/serovars.tsv with
    one row per .fasta in the -a folder. Set `write_results` to False to
    imitate a failed workflow, after which the real tool still exits 0.
    """
    from swineotype.adapters import app as app_mod

    tool = {"calls": [], "returncode": 0, "write_results": True}
    monkeypatch.setattr(app_mod.shutil, "which", lambda name: f"/opt/bin/{name}")

    def fake_run(cmd, cwd=None, **kwargs):
        staging = Path(cmd[cmd.index("-a") + 1])
        out = Path(cmd[cmd.index("-o") + 1])
        links = sorted(staging.iterdir())
        tool["calls"].append({"cmd": cmd, "cwd": Path(cwd), "staged": [p.name for p in links],
                              "targets": [p.resolve() for p in links]})
        if tool["write_results"]:
            write_serovars_tsv(out / "serovars.tsv", [(p.stem, "S1") for p in links])
        return subprocess.CompletedProcess(cmd, tool["returncode"])

    monkeypatch.setattr(app_mod.subprocess, "run", fake_run)
    return tool


def run_app(patterns, out_dir, summary=None):
    from swineotype.adapters import app as app_mod
    app_mod.run_app_analysis(assembly=[str(p) for p in patterns], out_dir=str(out_dir),
                             threads=2, swineotype_summary=str(summary) if summary else None)


def fasta_files(folder: Path, *names):
    folder.mkdir(parents=True, exist_ok=True)
    for name in names:
        (folder / name).write_text(">c\nACGT\n")
    return folder


def test_every_input_is_staged_as_sample_dot_fasta_and_handed_over(tmp_path, serovar_detector):
    inputs = fasta_files(tmp_path / "in", "iso1.fa", "iso2.fna", "iso3.fasta")

    run_app([inputs / "*"], tmp_path / "out")

    call, = serovar_detector["calls"]
    app_dir = (tmp_path / "out" / "app_detector").resolve()
    assert call["staged"] == ["iso1.fasta", "iso2.fasta", "iso3.fasta"]
    assert call["targets"] == [(inputs / n).resolve() for n in ("iso1.fa", "iso2.fna", "iso3.fasta")]
    assert call["cmd"] == ["/opt/bin/serovar_detector", "-a", str(app_dir / "assemblies"),
                           "-o", str(app_dir), "-t", "2"]
    assert call["cwd"] == app_dir, "Snakemake's .snakemake/ must land inside app_detector/"


def test_two_assembly_fasta_inputs_are_staged_under_distinct_names(tmp_path, serovar_detector):
    for run_name in ("runA", "runB"):
        fasta_files(tmp_path / "in" / run_name, "assembly.fasta")

    run_app([tmp_path / "in" / "*" / "assembly.fasta"], tmp_path / "out")

    staged = serovar_detector["calls"][0]["staged"]
    assert len(set(staged)) == 2, "two different files must not share one staged name"


def test_overlapping_patterns_stage_one_file_once(tmp_path, serovar_detector):
    inputs = fasta_files(tmp_path / "in", "iso1.fasta")

    run_app([inputs / "*.fasta", inputs / "iso1.fasta"], tmp_path / "out")

    assert serovar_detector["calls"][0]["staged"] == ["iso1.fasta"], "one file, one sample"


def test_a_rerun_hands_over_only_the_current_inputs(tmp_path, serovar_detector):
    """serovar_detector types every .fasta in the folder it is given, so a
    previous run's inputs left in it would be typed again."""
    inputs = fasta_files(tmp_path / "in", "iso1.fasta", "iso2.fasta", "iso3.fasta")

    run_app([inputs / "iso1.fasta", inputs / "iso2.fasta"], tmp_path / "out")
    run_app([inputs / "iso3.fasta"], tmp_path / "out")

    assert serovar_detector["calls"][1]["staged"] == ["iso3.fasta"]


def test_a_pattern_that_matched_nothing_is_rejected(tmp_path, serovar_detector):
    inputs = fasta_files(tmp_path / "in", "iso1.fasta")

    with pytest.raises(SystemExit) as exc:
        run_app([inputs / "*.fasta", inputs / "*.fna"], tmp_path / "out")
    assert exc.value.code == 2
    assert serovar_detector["calls"] == []


def test_a_failed_workflow_is_not_reported_as_success(tmp_path, serovar_detector):
    """serovar_detector exits 0 when its Snakemake workflow fails, so the exit
    code alone would report a run that produced nothing as a success."""
    inputs = fasta_files(tmp_path / "in", "iso1.fasta")
    serovar_detector["write_results"] = False

    with pytest.raises(SystemExit) as exc:
        run_app([inputs / "iso1.fasta"], tmp_path / "out")
    assert exc.value.code == 1


def test_a_previous_runs_table_does_not_pass_for_this_run(tmp_path, serovar_detector):
    inputs = fasta_files(tmp_path / "in", "iso1.fasta")
    run_app([inputs / "iso1.fasta"], tmp_path / "out")
    serovar_detector["write_results"] = False

    with pytest.raises(SystemExit) as exc:
        run_app([inputs / "iso1.fasta"], tmp_path / "out")
    assert exc.value.code == 1


def test_a_missing_serovar_detector_names_the_install_command(tmp_path, monkeypatch, capsys):
    from swineotype.adapters import app as app_mod
    monkeypatch.setattr(app_mod.shutil, "which", lambda name: None)
    inputs = fasta_files(tmp_path / "in", "iso1.fasta")

    with pytest.raises(SystemExit) as exc:
        run_app([inputs / "iso1.fasta"], tmp_path / "out")
    assert exc.value.code == 1
    assert "pip install ./third_party/serovar_detector" in capsys.readouterr().err


@pytest.mark.parametrize("out_name, file_name", [("my results", "iso1.fasta"),
                                                 ("out", "iso 1.fasta")])
def test_a_path_with_spaces_is_rejected_before_running(tmp_path, serovar_detector,
                                                        out_name, file_name):
    """serovar_detector matches paths with \\S+ and crashes on a space."""
    inputs = fasta_files(tmp_path / "in", file_name)

    with pytest.raises(SystemExit) as exc:
        run_app([inputs / file_name], tmp_path / out_name)
    assert exc.value.code == 2
    assert serovar_detector["calls"] == []


def test_the_calls_are_merged_into_the_summary(tmp_path, serovar_detector):
    summary = tmp_path / "all_isolates.csv"
    write_summary(summary, [{"sample": "suis1", "final_serotype": "2"}])
    inputs = fasta_files(tmp_path / "in", "app1.fasta")

    run_app([inputs / "app1.fasta"], tmp_path / "out", summary=summary)

    assert {r["sample"]: (r["final_serotype"], r["app_serovar"]) for r in read_rows(summary)} == \
        {"suis1": ("2", ""), "app1": ("", "S1")}


def test_the_run_credits_serovar_detector_and_records_both_versions(tmp_path, serovar_detector,
                                                                    monkeypatch, capsys):
    from swineotype import __version__
    from swineotype.adapters import app as app_mod
    monkeypatch.setattr(app_mod, "version", lambda dist: "1.1.2")
    inputs = fasta_files(tmp_path / "in", "app1.fasta")

    run_app([inputs / "app1.fasta"], tmp_path / "out")

    out = capsys.readouterr().out
    assert "serovar_detector 1.1.2 (Kasper Thystrup Karstensen)" in out
    assert "doi:10.1099/mgen.0.001434" in out
    assert "[OK] app1 => S1 (serovar_detector)" in out
    record = json.loads((tmp_path / "out" / "app_detector" / "swineotype_run.json").read_text())
    assert record["swineotype_version"] == __version__
    assert record["serovar_detector_version"] == "1.1.2"
    assert record["assemblies"] == [str((inputs / "app1.fasta").resolve())]
