
import os
import shutil
import tempfile
import pytest
from pathlib import Path
from swineotype.utils import ensure_unix_line_endings

def test_ensure_unix_line_endings():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        
        # Case 1: CRLF file
        crlf_file = tmp_path / "test_crlf.fasta"
        with open(crlf_file, "wb") as f:
            f.write(b">seq1\r\nATGC\r\n")
            
        fixed_path = ensure_unix_line_endings(str(crlf_file), str(tmp_path))
        
        assert fixed_path != str(crlf_file)
        assert Path(fixed_path).read_bytes() == b">seq1\nATGC\n"
        
        # Case 2: LF file
        lf_file = tmp_path / "test_lf.fasta"
        with open(lf_file, "wb") as f:
            f.write(b">seq1\nATGC\n")
            
        fixed_path_lf = ensure_unix_line_endings(str(lf_file), str(tmp_path))
        
        assert fixed_path_lf == str(lf_file)
        
if __name__ == "__main__":
    test_ensure_unix_line_endings()
    print("Test passed!")
