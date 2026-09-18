import hashlib
import os
import sys
from pathlib import Path

import click
import gzip

def ensure_tool(name: str):
    from shutil import which
    if which(name) is None:
        click.echo(f"[ERROR] Required tool not found in PATH: {name}", err=True)
        sys.exit(1)

def gzip_file(file_path: str):
    with open(file_path, 'rb') as f_in:
        with gzip.open(f"{file_path}.gz", 'wb') as f_out:
            f_out.writelines(f_in)
    os.remove(file_path)

def ensure_unix_line_endings(file_path: str, tmp_dir: str) -> str:
    """Stage an assembly into tmp_dir with LF line endings.

    We always stage a copy rather than reading in place, so that a read-only
    input mount (e.g. WSL) can still be processed.

    The staged name carries a short content digest. Staging on the bare
    basename meant two different assemblies called ``assembly.fasta`` clobbered
    each other in the shared tmp dir -- and, because the BLAST database and the
    samtools .fai index were also keyed on that name, one sample's coordinates
    could be applied to another sample's sequence.
    """
    path = Path(file_path)
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    content = path.read_bytes()
    # Normalize CRLF to LF, and also bare CR to LF just in case
    content = content.replace(b'\r\n', b'\n').replace(b'\r', b'\n')

    digest = hashlib.sha1(content).hexdigest()[:12]
    dest = tmp_dir / f"{path.stem}__{digest}{path.suffix}"

    # Same content already staged under this name: reuse it.
    if not (dest.exists() and dest.stat().st_size == len(content)):
        dest.write_bytes(content)

    return str(dest)
