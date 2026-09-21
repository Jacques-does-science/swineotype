"""CLI input handling and default persistence."""
import csv
import json
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from swineotype.main import DEFAULT_RESULTS_JSON, DEFAULT_RUN_JSON, DEFAULT_SUMMARY_CSV
from swineotype.main import SUMMARY_COLUMNS, expand_globs, main


def a_row(sample="test", status="STAGE1", final="1"):
    row = {c: "" for c in SUMMARY_COLUMNS}
    row.update({"sample": sample, "status": status, "final_serotype": final})
    return row


# --- glob expansion ----------------------------------------------------

def test_expand_globs_reports_a_pattern_that_matched_nothing(tmp_path):
    matched, unmatched = expand_globs([str(tmp_path / "*.fasta")])
    assert matched == []
    assert unmatched == [str(tmp_path / "*.fasta")]


def test_expand_globs_reports_a_missing_literal_path(tmp_path):
    matched, unmatched = expand_globs([str(tmp_path / "nope.fasta")])
    assert unmatched == [str(tmp_path / "nope.fasta")]


def test_expand_globs_reports_only_the_unmatched_pattern(tmp_path):
    (tmp_path / "a.fasta").write_text(">a\nACGT\n")
    matched, unmatched = expand_globs([str(tmp_path / "*.fasta"),
                                       str(tmp_path / "*.fna")])
    assert matched == [str(tmp_path / "a.fasta")]
    assert unmatched == [str(tmp_path / "*.fna")]


@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_an_unmatched_glob_exits_nonzero_without_writing_a_summary(mock_tool, mock_process):
    """Regression: an unmatched pattern expanded to nothing and the run
    finished successfully with an empty summary, so a typo'd path was
    indistinguishable from a batch where nothing was informative."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(main, ["--out_dir", "out", "--assembly", "missing_*.fasta"])
        assert result.exit_code != 0
        assert "No assembly matched" in result.output
        assert not Path("out", DEFAULT_SUMMARY_CSV).exists()
        mock_process.assert_not_called()


# --- default persistence -----------------------------------------------

@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_results_and_run_record_are_written_without_merged_csv(mock_tool, mock_process):
    """Regression: without --merged_csv the run left nothing machine-readable
    and nothing recording which thresholds and reference data produced it."""
    mock_process.return_value = a_row()
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("test.fasta").write_text(">t\nACGT\n")
        result = runner.invoke(main, ["--out_dir", "out", "--assembly", "test.fasta"])
        assert result.exit_code == 0, result.output

        summary = Path("out", DEFAULT_SUMMARY_CSV)
        assert summary.exists()
        rows = list(csv.DictReader(summary.open()))
        assert [r["sample"] for r in rows] == ["test"]
        assert list(rows[0]) == SUMMARY_COLUMNS

        results = json.loads(Path("out", DEFAULT_RESULTS_JSON).read_text())
        assert results[0]["final_serotype"] == "1"

        record = json.loads(Path("out", DEFAULT_RUN_JSON).read_text())
        assert record["swineotype_version"]
        assert record["status_counts"] == {"STAGE1": 1}
        assert record["effective_config"]["min_pid"] == 85.0
        assert record["effective_config"]["plurality"] == 0.60
        assert len(record["reference_data"]["resolver_refs_sha256"]) == 64
        assert record["reference_data"]["resolver_refs_manifest_version"]


@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_the_per_run_summary_is_rewritten_not_appended(mock_tool, mock_process):
    mock_process.return_value = a_row()
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("test.fasta").write_text(">t\nACGT\n")
        for _ in range(2):
            runner.invoke(main, ["--out_dir", "out", "--assembly", "test.fasta"])
        rows = list(csv.DictReader(Path("out", DEFAULT_SUMMARY_CSV).open()))
        assert len(rows) == 1, "the run's own summary reflects this run only"


@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_merged_csv_still_appends(mock_tool, mock_process):
    mock_process.return_value = a_row()
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("test.fasta").write_text(">t\nACGT\n")
        for _ in range(2):
            runner.invoke(main, ["--out_dir", "out", "--assembly", "test.fasta",
                                 "--merged_csv", "all.csv"])
        rows = list(csv.DictReader(Path("all.csv").open()))
        assert len(rows) == 2, "--merged_csv accumulates across runs, as documented"


@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_input_species_is_passed_through(mock_tool, mock_process):
    mock_process.return_value = a_row()
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("test.fasta").write_text(">t\nACGT\n")
        runner.invoke(main, ["--out_dir", "out", "--assembly", "test.fasta",
                             "--input_species", "Streptococcus suis"])
        assert mock_process.call_args.kwargs["input_species"] == "Streptococcus suis"
        record = json.loads(Path("out", DEFAULT_RUN_JSON).read_text())
        assert record["input_species"] == "Streptococcus suis"


@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_appending_to_a_summary_with_other_columns_is_refused(mock_tool, mock_process):
    """Regression: DictWriter happily writes new-shaped rows under an old
    header, putting every value in the wrong column. The summary gained and
    lost columns in 0.2.0, so an older merged CSV must be refused, loudly."""
    mock_process.return_value = a_row()
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("test.fasta").write_text(">t\nACGT\n")
        # A merged CSV in the pre-0.2.0 shape.
        Path("old.csv").write_text(
            "sample,sample_path,run_dir,species,stage1_top,status,final_serotype\n"
            "iso1,/x/iso1.fasta,iso1,Streptococcus suis,2,STAGE2,2\n")

        result = runner.invoke(main, ["--out_dir", "out", "--assembly", "test.fasta",
                                      "--merged_csv", "old.csv"])

        assert result.exit_code != 0
        assert "different column set" in result.output
        assert Path("old.csv").read_text().count("\n") == 2, "the old file is untouched"
        # The run's own outputs are still written: only the append is refused.
        assert Path("out", DEFAULT_SUMMARY_CSV).exists()


@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_appending_to_a_matching_summary_is_allowed(mock_tool, mock_process):
    mock_process.return_value = a_row()
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("test.fasta").write_text(">t\nACGT\n")
        for _ in range(3):
            r = runner.invoke(main, ["--out_dir", "out", "--assembly", "test.fasta",
                                     "--merged_csv", "all.csv"])
            assert r.exit_code == 0, r.output
        assert len(list(csv.DictReader(Path("all.csv").open()))) == 3
