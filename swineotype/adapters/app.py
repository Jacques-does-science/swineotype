#!/usr/bin/env python3
"""APP serotyping, by serovar_detector.

serovar_detector (Kasper Thystrup Karstensen, MIT licence) does the APP typing.
It is bundled unmodified as a git submodule pinned to one of its releases and is
installed from there. swineotype only stages the assemblies, runs
serovar_detector's own command-line tool on them, and joins the serovar it
suggests onto the S. suis summary.
"""

import csv
import datetime as _dt
import json
import re
import shutil
import subprocess
import sys
from glob import glob
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Optional, List

import click

from swineotype import __version__
from swineotype.utils import unique_run_names

CITATION = ("Angen Ø, Karstensen KT, Vilaró A, et al. (2025) Serotyping of Actinobacillus "
            "pleuropneumoniae based on whole genome sequencing: validation of a bioinformatic "
            "tool. Microb Genom 11(7):001434. doi:10.1099/mgen.0.001434")


def log(msg: str):
    click.echo(f"[INFO] {msg}")

def err(msg: str):
    click.echo(f"[ERROR] {msg}", file=sys.stderr)


def read_swineotype_summary(path: Path) -> tuple[list[str], list[dict]]:
    """Read the CSV summary `swineotype --species suis` writes: (columns, rows).

    It used to be read with ``sep="\\t"``, which does not raise on a CSV:
    pandas returned one column named after the whole header line, so the
    try/except around it never fired and the merge died with
    ``KeyError("sample")``.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        columns = reader.fieldnames  # read while the file is open: it is lazy
        rows = list(reader)
    if not columns:
        raise ValueError(f"{path} is empty; expected a summary written by `swineotype --species suis`")
    if "sample" not in columns:
        raise ValueError(f"{path} has no 'sample' column (found: {columns}); "
                         f"expected a summary written by `swineotype --species suis`")
    return list(columns), rows


def read_app_calls(app_results: Path) -> dict[str, str]:
    """Sample -> suggested serovar, from serovar_detector's serovars.tsv."""
    with open(app_results, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"Sample", "Suggested_serovar"} - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{app_results} is missing expected column(s): {sorted(missing)}")
        return {r["Sample"]: r["Suggested_serovar"] for r in reader}


def merge_app_results(summary: Path, app_results: Path) -> None:
    """Outer-join the APP calls onto the swineotype summary by `sample`, in place.

    A merge is ordered by sample and a fresh summary keeps serovars.tsv's order.
    Values are copied as text: pandas re-typed any all-numeric column, so an
    outer join that introduced a blank turned `final_serotype` "2" into "2.0".
    """
    calls = read_app_calls(app_results)
    if summary.exists():
        columns, rows = read_swineotype_summary(summary)
        seen = {r["sample"] for r in rows}
        rows = [{**r, "app_serovar": calls.get(r["sample"], r.get("app_serovar", ""))} for r in rows]
        rows = sorted(rows + [{"sample": s, "app_serovar": v} for s, v in calls.items() if s not in seen],
                      key=lambda r: r["sample"])
    else:
        columns, rows = ["sample"], [{"sample": s, "app_serovar": v} for s, v in calls.items()]
    if "app_serovar" not in columns:
        columns = columns + ["app_serovar"]
    with open(summary, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run_app_analysis(assembly: List[str], out_dir: str, threads: int, swineotype_summary: Optional[str]):
    """Stage the assemblies, run serovar_detector on them, and merge its calls."""
    app_dir = Path(out_dir).resolve() / "app_detector"

    # Expand absolute glob patterns safely, and fail on any pattern that
    # matched nothing rather than quietly processing the ones that did.
    assemblies = []
    unmatched = []
    for p in assembly:
        pattern = p.strip('"').strip("'")
        matches = sorted(glob(pattern))
        if not matches:
            unmatched.append(pattern)
        assemblies.extend(Path(g).resolve() for g in matches)
    # Overlapping patterns must not stage one file twice.
    assemblies = list(dict.fromkeys(assemblies))

    if unmatched:
        for pattern in unmatched:
            err(f"No assembly matched: {pattern}")
        sys.exit(2)
    if not assemblies:
        err("No assemblies to process")
        sys.exit(2)
    log(f"Found {len(assemblies)} assemblies")

    # Each file is linked as `{sample}.fasta` in one folder: serovar_detector
    # names a sample after its file, and its -a mode types every .fasta/.fa in
    # the folder it is given.
    staging = app_dir / "assemblies"
    staged = {staging / f"{name}.fasta": asm_path
              for asm_path, name in zip(assemblies, unique_run_names(assemblies))}
    # It matches file paths with \S+ and puts them into shell commands
    # unquoted, so a path containing whitespace makes it crash.
    spaced = [str(link) for link in staged if re.search(r"\s", str(link))]
    if spaced:
        for link in spaced:
            err(f"serovar_detector cannot handle paths that contain spaces: {link}")
        err("Rename those files, or choose an --out_dir without spaces.")
        sys.exit(2)

    executable = shutil.which("serovar_detector")
    if executable is None:
        err("serovar_detector is not installed. From the swineotype checkout, install the "
            "bundled release with `pip install ./third_party/serovar_detector` "
            "(scripts/install_swineotype.sh does this).")
        sys.exit(1)
    try:
        sd_version = version("serovar_detector")
    except PackageNotFoundError:
        sd_version = "unknown"
    log(f"APP serovars are called by serovar_detector {sd_version} "
        f"(Kasper Thystrup Karstensen). If you use them, please cite: {CITATION}")

    # Rebuilt on every run, so that the folder holds exactly this run's inputs.
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True)
    for link, asm_path in staged.items():
        link.symlink_to(asm_path)

    # serovar_detector exits 0 even when its workflow fails, so success is the
    # results table appearing -- which a previous run's copy must not fake.
    app_results = app_dir / "serovars.tsv"
    app_results.unlink(missing_ok=True)

    cmd = [executable, "-a", str(staging), "-o", str(app_dir), "-t", str(threads)]
    (app_dir / "swineotype_run.json").write_text(json.dumps({
        "swineotype_version": __version__,
        "serovar_detector_version": sd_version,
        "run_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "command": cmd,
        "assemblies": [str(a) for a in assemblies],
    }, indent=2) + "\n")
    log(f"Command: {' '.join(cmd)}")
    # serovar_detector runs Snakemake in the current directory, which is where
    # Snakemake keeps its .snakemake/ bookkeeping: keep that inside app_dir.
    ret = subprocess.run(cmd, cwd=app_dir)
    if ret.returncode != 0 or not app_results.exists():
        err(f"serovar_detector did not produce {app_results}; see its output above.")
        sys.exit(ret.returncode or 1)

    for sample, serovar in read_app_calls(app_results).items():
        click.echo(f"[OK] {sample} => {serovar} (serovar_detector)")
    log(f"APP serovar table: {app_results}")

    # Optional merge with the swineotype summary
    if swineotype_summary:
        swineo = Path(swineotype_summary).resolve()
        log(f"Merging APP results into {swineo}")
        try:
            merge_app_results(swineo, app_results)
        except ValueError as exc:
            err(str(exc))
            sys.exit(1)
        log(f"[SUCCESS] Summary written → {swineo}")
