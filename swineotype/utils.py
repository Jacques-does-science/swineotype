import sys
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
    import os
    os.remove(file_path)

def ensure_unix_line_endings(file_path: str, tmp_dir: str) -> str:
    from pathlib import Path
    path = Path(file_path)
    dest = Path(tmp_dir) / path.name

    # Even if line endings are fine, we copy to tmp_dir 
    # to ensure we have write permission for the .fai index file.
    # This addresses issues where input is in a read-only mount (e.g. WSL).
    
    # If source and dest resolve to the same file (e.g. user provided file inside tmp_dir),
    # we might overwrite it. To be safe, use a prefixed name if they are the same.
    if dest.resolve() == path.resolve():
         dest = Path(tmp_dir) / f"staged_{path.name}"
    
    # click.echo(f"[INFO] Staging assembly to {dest}...")
    
    with open(path, 'rb') as f_in:
        content = f_in.read()
        
    with open(dest, 'wb') as f_out:
        # Normalize CRLF to LF, and also bare CR to LF just in case
        content = content.replace(b'\r\n', b'\n').replace(b'\r', b'\n')
        f_out.write(content)
        
    return str(dest)
