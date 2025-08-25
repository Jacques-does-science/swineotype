from unittest.mock import patch
from swineotype.blast import run_blast

@patch("subprocess.run")
def test_run_blast(mock_run):
    run_blast("query.fasta", "db_prefix", 4, "outfmt")
    mock_run.assert_called_once_with(
        [
            "blastn",
            "-query",
            "query.fasta",
            "-db",
            "db_prefix",
            "-task",
            "blastn",
            "-outfmt",
            "outfmt",
            "-max_target_seqs",
            "50",
            "-num_threads",
            "4",
        ],
        check=True,
        capture_output=True,
        cwd=None,
        text=True,
    )
