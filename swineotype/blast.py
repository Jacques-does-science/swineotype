import shlex
import subprocess
from pathlib import Path

def run(cmd, check=True, capture=True, cwd=None, text=True):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    res = subprocess.run(cmd, check=check, capture_output=capture, cwd=cwd, text=text)
    return res.stdout

def make_db_if_needed(asm_fa: str, tmpdir: Path) -> str:
    prefix = tmpdir / ("asmdb_" + Path(asm_fa).stem)
    nin = prefix.with_suffix(".nin")
    ndb = prefix.with_suffix(".ndb")
    if not (nin.exists() or ndb.exists()):
        run(["makeblastdb", "-in", asm_fa, "-dbtype", "nucl", "-out", str(prefix)])
    return str(prefix)

def run_blast(query_fa: str, db_prefix: str, threads: int, outfmt_cols: str, max_target_seqs=50) -> str:
    cmd = [
        "blastn",
        "-query", query_fa,
        "-db", db_prefix,
        "-task", "blastn",
        "-outfmt", outfmt_cols,
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
    ]
    return run(cmd)
