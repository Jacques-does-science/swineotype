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
        result = runner.invoke(main, ["--out_dir", "out", "test.fasta"])
        assert result.exit_code == 0
        assert "[OK] test.fasta => 1 (STAGE1)" in result.output
        mock_ensure_tool.assert_any_call("blastn")
        mock_process_one.assert_called_once()
