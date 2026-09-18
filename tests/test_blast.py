from pathlib import Path
from unittest.mock import patch

from swineotype.blast import db_prefix_for, run_blast


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
            "-dust",
            "no",
        ],
        check=True,
        capture_output=True,
        cwd=None,
        text=True,
    )


def test_db_prefix_distinguishes_same_basename(tmp_path):
    """Regression: two isolates both called assembly.fasta must not share a DB.

    Flye writes assembly.fasta and SPAdes writes contigs.fasta, so batch runs
    over per-sample output directories hit this constantly. Keying the cache on
    the basename alone meant the second sample was BLASTed against the first
    sample's database.
    """
    a = tmp_path / "runA" / "assembly.fasta"
    b = tmp_path / "runB" / "assembly.fasta"
    a.parent.mkdir(); b.parent.mkdir()
    a.write_text(">c1\nACGTACGTAA\n")
    b.write_text(">c1\nTTGCATGCAT\n")

    cache = tmp_path / "cache"; cache.mkdir()
    assert db_prefix_for(str(a), cache) != db_prefix_for(str(b), cache)


def test_db_prefix_is_stable_for_identical_content(tmp_path):
    """Same content under different names may safely share a cached DB."""
    a = tmp_path / "one.fasta"
    b = tmp_path / "one.fasta.copy"
    a.write_text(">c1\nACGT\n")
    b.write_text(">c1\nACGT\n")
    cache = tmp_path / "cache"; cache.mkdir()
    # differing stems, identical digests
    assert db_prefix_for(str(a), cache).split("_")[-1] == db_prefix_for(str(b), cache).split("_")[-1]


def test_db_prefix_survives_dots_in_stem(tmp_path):
    """Path.with_suffix() would have truncated 'asmdb_sample.v2' to 'asmdb_sample'."""
    f = tmp_path / "sample.v2.fasta"
    f.write_text(">c1\nACGT\n")
    cache = tmp_path / "cache"; cache.mkdir()
    prefix = db_prefix_for(str(f), cache)
    assert Path(prefix).name.startswith("asmdb_sample.v2_")
