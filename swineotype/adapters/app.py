#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path
from glob import glob
from typing import Optional, List

import pandas as pd
import yaml
import click

def log(msg: str):
    click.echo(f"[INFO] {msg}")

def err(msg: str):
    click.echo(f"[ERROR] {msg}", file=sys.stderr)


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

    # --- FIX: expand absolute glob patterns safely ---
    assemblies = []
    for p in assembly:
        pattern = p.strip('"').strip("'")
        assemblies.extend(Path(g).resolve() for g in glob(pattern))

    if not assemblies:
        err(f"No assemblies found for pattern(s): {', '.join(assembly)}")
        sys.exit(1)
    log(f"Found {len(assemblies)} assemblies")

    # Write sample_sheet.csv for Snakemake/peppy
    sample_sheet_csv = schemas_dir / "sample_sheet.csv"
    with open(sample_sheet_csv, "w") as fh:
        fh.write("sample_name,type\n")
        for fa in assemblies:
            fh.write(f"{fa.stem},Assembly\n")
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
    for asm_path in assemblies:
        symlink_path = assembly_tmp_dir / asm_path.name
        if not symlink_path.exists():
            symlink_path.symlink_to(asm_path)
        else:
            # Overwrite if it's a broken link, for example
            if not symlink_path.resolve(strict=False).exists():
                symlink_path.unlink()
                symlink_path.symlink_to(asm_path)


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
        
        app_df = pd.read_csv(app_results, sep="\t")
        app_clean = app_df.rename(
            columns={"Sample": "sample", "Suggested_serovar": "app_serovar"}
        )[["sample", "app_serovar"]]

        if swineo.exists():
            log(f"Merging APP results with existing summary → {swineo}")
            # Try TSV then CSV automatically
            try:
                suis_df = pd.read_csv(swineo, sep="\t")
            except Exception:
                suis_df = pd.read_csv(swineo)

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
