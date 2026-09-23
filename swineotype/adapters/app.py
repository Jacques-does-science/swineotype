#!/usr/bin/env python3

import csv
import subprocess
import sys
from pathlib import Path
from glob import glob
from typing import Optional, List

import yaml
import click

from swineotype import __version__
from swineotype.utils import unique_run_names


def log(msg: str):
    click.echo(f"[INFO] {msg}")

def err(msg: str):
    click.echo(f"[ERROR] {msg}", file=sys.stderr)


def stage_assembly(source: Path, dest: Path) -> None:
    """Point `dest` at `source`, replacing whatever is already there.

    `Path.exists()` follows symlinks, so a dangling link answered False and
    the `symlink_to()` that followed raised FileExistsError -- while a link
    that still resolved was kept even when it pointed at a previous run's file
    of the same name. Both are fixed by never reusing an existing entry.
    """
    dest.unlink(missing_ok=True)
    dest.symlink_to(source)


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


def merge_app_results(summary: Path, app_results: Path) -> None:
    """Outer-join the APP calls onto the swineotype summary by `sample`, in place.

    serovar.tsv is written by readr::write_tsv in the workflow's R summariser:
    a real TSV with a known header. A merge is ordered by sample and a fresh
    summary keeps serovar.tsv's order, as the pandas code this replaces did.
    Values are copied as text: pandas re-typed any all-numeric column, so an
    outer join that introduced a blank turned `final_serotype` "2" into "2.0".
    """
    with open(app_results, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"Sample", "Suggested_serovar"} - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{app_results} is missing expected column(s): {sorted(missing)}")
        calls = {r["Sample"]: r["Suggested_serovar"] for r in reader}

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

    """Adapter for APP serovar detection + merge with swineotype"""
    outdir = Path(out_dir).resolve()
    app_dir = outdir / "app_detector"
    results_dir = app_dir / "results"
    tmp_dir = app_dir / "tmp"
    config_dir = app_dir / "config"
    logs_dir = app_dir / "logs"
    schemas_dir = results_dir / "schemas"

    for d in (results_dir, tmp_dir, config_dir, logs_dir, schemas_dir):
        d.mkdir(parents=True, exist_ok=True)

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
    # Overlapping patterns must not produce two sample-sheet rows and two
    # staged links for one file.
    assemblies = list(dict.fromkeys(assemblies))

    if unmatched:
        for pattern in unmatched:
            err(f"No assembly matched: {pattern}")
        sys.exit(2)
    if not assemblies:
        err("No assemblies to process")
        sys.exit(2)
    log(f"Found {len(assemblies)} assemblies")

    sample_names = unique_run_names(assemblies)

    # Write sample_sheet.csv for Snakemake/peppy
    sample_sheet_csv = schemas_dir / "sample_sheet.csv"
    with open(sample_sheet_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_name", "type"])
        for name in sample_names:
            writer.writerow([name, "Assembly"])
    log(f"Wrote samples table with {len(assemblies)} assemblies → {sample_sheet_csv}")


    # KMA DB prefix (must exist): .../third_party/serovar_detector/db/Actinobacillus_pleuropneumoniae.*
    third_party = Path(__file__).parent.parent.parent / "third_party" / "serovar_detector"
    db_dir = third_party / "db"
    db_prefix = db_dir / "Actinobacillus_pleuropneumoniae"
    if not all(db_prefix.with_suffix(s).exists() for s in (".fasta", ".seq.b", ".comp.b", ".length.b")):
        err(f"Database prefix not found or incomplete: {db_prefix}")
        sys.exit(1)
    log(f"Using KMA DB prefix: {db_prefix}")

    serovar_profiles = third_party / "config" / "serovar_profiles.yaml"
    if not serovar_profiles.exists():
        err(f"Missing serovar profiles YAML: {serovar_profiles}")
        sys.exit(1)

    # Copy serovar_profiles to config_dir for the R script
    import shutil
    shutil.copy(serovar_profiles, config_dir / "serovar_profiles.yaml")

    # Create the peppy project config, which will live in the schemas_dir
    project_cfg = schemas_dir / "project_config.yaml"
    if not project_cfg.exists():
        project_cfg.write_text(
            "pep_version: 2.1.0\n"
            "name: app_serovar_project\n"
            "sample_table: sample_sheet.csv\n"  # Points to the CSV in the same directory
        )

    config_yaml = config_dir / "config.yaml"
    config = {
        "outdir": str(results_dir),
        "tmpdir": str(tmp_dir),
        "append_results": False,
        "database": str(db_prefix),
        "threads": int(threads),
        "threshold": 98.0,
        "debug": False,
        "serovar_profiles": str(serovar_profiles),
        "summary_file": str(results_dir / "serovar_summary.tsv"),
        "log_dir": str(logs_dir),
        "results_dir": str(results_dir),
        "schemas": str(schemas_dir),
        "version": "2.1.0",
        "swineotype_version": __version__,
    }
    with open(config_yaml, "w") as fh:
        yaml.dump(config, fh)
    log(f"Workflow config written: {config_yaml}")

    # Create symlinks for the assembly files in the tmp directory, which is what
    # the serovar_detector workflow expects.
    assembly_tmp_dir = tmp_dir / "assemblies"
    assembly_tmp_dir.mkdir(parents=True, exist_ok=True)
    log(f"Creating symlinks for assemblies in {assembly_tmp_dir}")
    for asm_path, name in zip(assemblies, sample_names):
        # `{sample}.fasta`, always: the workflow's input pattern is literal, so
        # staging a .fa or .fna under its own name left the rule with no input.
        stage_assembly(asm_path, assembly_tmp_dir / f"{name}.fasta")


    # Run Snakemake
    snakefile = third_party / "workflow" / "Snakefile"
    cmd = [
        "snakemake",
        "-s", str(snakefile),
        "--configfile", str(config_yaml),
        "--cores", str(threads),
        "--directory", str(app_dir),
        "--use-conda",
    ]
    log(f"Command: {' '.join(cmd)}")
    ret = subprocess.run(cmd)
    if ret.returncode != 0:
        err("SerovarDetector failed.")
        sys.exit(ret.returncode)

    # APP results
    app_results = results_dir / "serovar.tsv"
    if not app_results.exists():
        err(f"APP serovar results not found: {app_results}")
        sys.exit(1)

    # Optional merge with the swineotype summary
    if swineotype_summary:
        swineo = Path(swineotype_summary).resolve()
        log(f"Merging APP results into {swineo}")
        try:
            merge_app_results(swineo, app_results)
        except ValueError as exc:
            err(str(exc))
            sys.exit(1)
        log(f"[SUCCESS] Summary written \u2192 {swineo}")
