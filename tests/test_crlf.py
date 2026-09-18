import tempfile
from pathlib import Path

from swineotype.utils import ensure_unix_line_endings


def test_crlf_is_normalised():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        crlf_file = tmp_path / "test_crlf.fasta"
        crlf_file.write_bytes(b">seq1\r\nATGC\r\n")

        staged = ensure_unix_line_endings(str(crlf_file), str(tmp_path))

        assert staged != str(crlf_file)
        assert Path(staged).read_bytes() == b">seq1\nATGC\n"


def test_lf_input_is_still_staged():
    """Input is always copied, even when already LF, so a read-only input
    mount (e.g. WSL) can still be processed."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        lf_file = tmp_path / "test_lf.fasta"
        lf_file.write_bytes(b">seq1\nATGC\n")

        staged = ensure_unix_line_endings(str(lf_file), str(tmp_path))

        assert staged != str(lf_file)
        assert Path(staged).read_bytes() == b">seq1\nATGC\n"


def test_same_basename_different_content_does_not_clobber(tmp_path):
    """Regression: staging on the bare basename let two different assemblies
    named assembly.fasta overwrite each other in the shared tmp dir."""
    a = tmp_path / "runA" / "assembly.fasta"
    b = tmp_path / "runB" / "assembly.fasta"
    a.parent.mkdir(); b.parent.mkdir()
    a.write_bytes(b">c1\nAAAA\n")
    b.write_bytes(b">c1\nCCCC\n")
    staging = tmp_path / "staging"

    staged_a = ensure_unix_line_endings(str(a), str(staging))
    staged_b = ensure_unix_line_endings(str(b), str(staging))

    assert staged_a != staged_b
    assert Path(staged_a).read_bytes() == b">c1\nAAAA\n"
    assert Path(staged_b).read_bytes() == b">c1\nCCCC\n"
