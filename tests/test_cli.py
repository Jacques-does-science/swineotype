from unittest.mock import patch
from click.testing import CliRunner
from swineotype.main import main

@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_main_cli(mock_ensure_tool, mock_process_one):
    mock_process_one.return_value = {
        "sample": "test.fasta",
        "status": "STAGE1",
        "final_serotype": "1",
    }
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("test.fasta", "w") as f:
            f.write(">test\nACGT")
        result = runner.invoke(main, ["--out_dir", "out", "--assembly", "test.fasta"])
        assert result.exit_code == 0
        assert "[OK] test.fasta => 1 (STAGE1)" in result.output
        mock_ensure_tool.assert_any_call("blastn")
        mock_process_one.assert_called_once()


@patch("swineotype.main.process_one")
@patch("swineotype.main.ensure_tool")
def test_main_cli_multiple_assemblies(mock_ensure_tool, mock_process_one):
    mock_process_one.return_value = {
        "sample": "test.fasta",
        "status": "STAGE1",
        "final_serotype": "1",
    }
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("a.fasta", "w") as f:
            f.write(">a\nACGT")
        with open("b.fasta", "w") as f:
            f.write(">b\nACGT")
        result = runner.invoke(main, ["--out_dir", "out", "--assembly", "a.fasta", "--assembly", "b.fasta"])
        assert result.exit_code == 0
        assert mock_process_one.call_count == 2


@patch("swineotype.main.run_app_analysis")
def test_main_cli_app(mock_run_app_analysis):
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("a.fasta", "w") as f:
            f.write(">a\nACGT")
        with open("b.fasta", "w") as f:
            f.write(">b\nACGT")
        result = runner.invoke(main, ["--species", "app", "--out_dir", "out", "--assembly", "*.fasta"])
        assert result.exit_code == 0
        mock_run_app_analysis.assert_called_once()
        # The first argument of the first call to the mock
        called_args, called_kwargs = mock_run_app_analysis.call_args
        assert called_kwargs['assembly'] == ['*.fasta']
