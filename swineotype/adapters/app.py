#!/usr/bin/env python3

import csv
import subprocess
import sys
from pathlib import Path
from glob import glob
from typing import Optional, List

import pandas as pd
import yaml
import click

from swineotype import __version__
from swineotype.utils import unique_run_names

# The workflow addresses every staged assembly as {sample}.fasta, so the
# staged copy has to carry that exact suffix whatever the input was called.
STAGED_SUFFIX = ".fasta"


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


def read_swineotype_summary(path: Path) -> pd.DataFrame:
    """Read the CSV summary `swineotype --species suis` writes.

    It used to be read with ``sep="\\t"``, which does not raise on a CSV:
    pandas returns one column named after the whole header line, so the
    try/except around it never fired and the merge died with
    ``KeyError("sample")``. An empty file raises EmptyDataError, a ValueError.
    """
    df = pd.read_csv(path, dtype={"sample": str})
    if "sample" not in df.columns:
        raise ValueError(f"{path} has no 'sample' column (found: {list(df.columns)}); "
                         f"expected a summary written by `swineotype --species suis`")
    return df


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
    if not (db_prefix.with_suffix(".fasta").exists()
            and db_prefix.with_suffix(".seq.b").exists()
            and db_prefix.with_suffix(".comp.b").exists()
            and db_prefix.with_suffix(".length.b").exists()):
        err(f"Database prefix not found or incomplete: {db_prefix}")
        sys.exit(1)
    log(f"Using KMA DB prefix: {db_prefix}")

    # Paths in config
    # ...
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

    log("[INFO] Effective config.yaml contents:")
    print(yaml.dump(config, sort_keys=False))

    # Create symlinks for the assembly files in the tmp directory, which is what
    # the serovar_detector workflow expects.
    assembly_tmp_dir = tmp_dir / "assemblies"
    assembly_tmp_dir.mkdir(parents=True, exist_ok=True)
    log(f"Creating symlinks for assemblies in {assembly_tmp_dir}")
    for asm_path, name in zip(assemblies, sample_names):
        # `{sample}.fasta`, always: the workflow's input pattern is literal, so
        # staging a .fa or .fna under its own name left the rule with no input.
        stage_assembly(asm_path, assembly_tmp_dir / f"{name}{STAGED_SUFFIX}")


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

    # Optional merge with swineotype summary
    # Optional merge with swineotype summary or just output to CSV
    if swineotype_summary:
        swineo = Path(swineotype_summary).resolve()
        
        # serovar.tsv is written by readr::write_tsv in the workflow's R
        # summariser: a real TSV, with a known header.
        app_df = pd.read_csv(app_results, sep="\t", dtype=str)
        missing = {"Sample", "Suggested_serovar"} - set(app_df.columns)
        if missing:
            err(f"{app_results} is missing expected column(s): {sorted(missing)}")
            sys.exit(1)
        app_clean = app_df.rename(
            columns={"Sample": "sample", "Suggested_serovar": "app_serovar"}
        )[["sample", "app_serovar"]]

        if swineo.exists():
            log(f"Merging APP results with existing summary → {swineo}")
            try:
                suis_df = read_swineotype_summary(swineo)
            except ValueError as exc:
                err(str(exc))
                sys.exit(1)

            merged = suis_df.merge(app_clean, on="sample", how="outer")
            merged.to_csv(swineo, index=False)
            log(f"[SUCCESS] Updated summary written → {swineo}")
        else:
            log(f"Writing APP results to new summary → {swineo}")
            app_clean.to_csv(swineo, index=False)
            log(f"[SUCCESS] Summary written → {swineo}")

@click.command()
@click.option("--assembly", multiple=True, required=True, help="Path to one or more assembly files or glob patterns.")
@click.option("--out_dir", required=True, help="Output directory base")
@click.option("--threads", type=int, default=4, help="Threads for Snakemake/KMA")
@click.option("--swineotype_summary", help="Path to swineotype summary TSV/CSV to merge with APP results")
def main(assembly, out_dir, threads, swineotype_summary):
    """Adapter for APP serovar detection + merge with swineotype"""
    run_app_analysis(
        assembly=list(assembly),
        out_dir=out_dir,
        threads=threads,
        swineotype_summary=swineotype_summary
    )

if __name__ == "__main__":
    main()
