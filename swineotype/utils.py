import sys
import click

def ensure_tool(name: str):
    from shutil import which
    if which(name) is None:
        click.echo(f"[ERROR] Required tool not found in PATH: {name}", err=True)
        sys.exit(1)
