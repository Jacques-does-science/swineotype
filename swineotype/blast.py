import hashlib
import subprocess
from pathlib import Path

def run(cmd: list[str]) -> str:
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout

def content_digest(path, length: int = 12) -> str:
    """Short SHA-1 of a file's bytes, used to key caches by content."""
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:length]

def db_prefix_for(asm_fa: str, tmpdir: Path) -> str:
    """Cache key for an assembly's BLAST database.

    Keyed on file CONTENT, not on the basename. Keying on
    ``Path(asm_fa).stem`` alone meant two different assemblies that happen to
    share a filename -- ``assembly.fasta`` from Flye, ``contigs.fasta`` from
    SPAdes, i.e. the common case when batch-processing per-sample output
    directories -- collided in the shared cache, and the second sample was
    silently BLASTed against the first sample's database.
    """
    return str(Path(tmpdir) / f"asmdb_{Path(asm_fa).stem}_{content_digest(asm_fa)}")

def make_db_if_needed(asm_fa: str, tmpdir: Path) -> str:
    prefix = db_prefix_for(asm_fa, tmpdir)
    # Build the sibling paths by string concatenation: Path.with_suffix() would
    # eat the last dot-separated segment of stems like "sample.v2".
    if not (Path(prefix + ".nin").exists() or Path(prefix + ".ndb").exists()):
        run(["makeblastdb", "-in", asm_fa, "-dbtype", "nucl", "-out", prefix])
    return prefix

def run_blast(query_fa: str, db_prefix: str, threads: int, outfmt_cols: str, max_target_seqs=50) -> str:
    cmd = [
        "blastn",
        "-query", query_fa,
        "-db", db_prefix,
        "-task", "blastn",
        "-outfmt", outfmt_cols,
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
        "-dust", "no",
    ]
    return run(cmd)
