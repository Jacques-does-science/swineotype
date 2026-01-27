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
    # Check if we have CRLF
    has_crlf = False
    with open(path, 'rb') as f:
        while chunk := f.read(4096):
            if b'\r' in chunk:
                has_crlf = True
                break
    
    if not has_crlf:
        return file_path
    
    # Needs conversion
    dest = Path(tmp_dir) / path.name
    
    # If source and dest are the same, rename dest to avoid overwriting/truncation issues
    if dest.resolve() == path.resolve():
         dest = Path(tmp_dir) / f"{path.stem}_unix{path.suffix}"

    # Avoid overwriting if unrelated file exists, but good to refresh
    click.echo(f"[INFO] Detected Windows line endings in {path.name}. Creating normalized copy at {dest}...")
    
    with open(path, 'rb') as f_in:
        content = f_in.read()
        
    with open(dest, 'wb') as f_out:
        content = content.replace(b'\r\n', b'\n')
        f_out.write(content)
        
    return str(dest)
