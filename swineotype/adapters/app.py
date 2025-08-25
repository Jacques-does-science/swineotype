#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path
from glob import glob
import pandas as pd
import yaml
import click

def log(msg: str):
    click.echo(f"[INFO] {msg}", flush=True)

def err(msg: str):
    click.echo(f"[ERROR] {msg}", file=sys.stderr, flush=True)

@click.command()
@click.argument("assembly", required=True)
@click.option("--out_dir", required=True, help="Output directory base")
@click.option("--threads", type=int, default=4, help="Threads for Snakemake/KMA")
@click.option("--swineotype_summary", help="Path to swineotype summary TSV/CSV to merge with APP results")
def main(assembly, out_dir, threads, swineotype_summary):
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
    pattern = assembly.strip('"').strip("'")
    assemblies = [Path(p).resolve() for p in glob(pattern)]
    if not assemblies:
        err(f"No assemblies found for pattern: {assembly}")
        sys.exit(1)
    log(f"Found {len(assemblies)} assemblies")

    # Write samples.tsv for Snakemake
    samples_tsv = app_dir / "samples.tsv"
    with open(samples_tsv, "w") as fh:
        fh.write("sample_name\tassembly\ttype\n")
        for fa in assemblies:
            fh.write(f"{fa.stem}\t{fa}\tAssembly\n")
    log(f"Wrote samples table with {len(assemblies)} assemblies → {samples_tsv}")

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
    serovar_profiles = third_party / "config" / "serovar_profiles.yaml"
    if not serovar_profiles.exists():
        err(f"Missing serovar profiles YAML: {serovar_profiles}")
        sys.exit(1)

    project_cfg = schemas_dir / "project_config.yaml"
    if not project_cfg.exists():
        # Minimal PEP config to satisfy Snakefile/peppy, not actually used for samples here.
        project_cfg.write_text(
            "pep_version: 2.1.0\n"
            "name: app_serovar_project\n"
            "samples: ../samples.tsv\n"
        )

    config_yaml = config_dir / "config.yaml"
    config = {
        "outdir": str(results_dir),
        "tmpdir": str(tmp_dir),
        "append_results": False,
        "database": str(db_prefix),
        "samples": str(samples_tsv),
        "threads": int(threads),
        "threshold": 98.0,
        "debug": False,
        "project_config": str(project_cfg),
        "type": "Assembly",
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
    if swineotype_summary:
        swineo = Path(swineotype_summary).resolve()
        if not swineo.exists():
            err(f"Swineotype summary not found: {swineo}")
            sys.exit(1)

        log(f"Merging APP results with swineotype summary → {swineo}")
        # Try TSV then CSV automatically
        try:
            suis_df = pd.read_csv(swineo, sep="\t")
        except Exception:
            suis_df = pd.read_csv(swineo)

        app_df = pd.read_csv(app_results, sep="\t")

        app_clean = app_df.rename(
            columns={"Sample": "sample", "Suggested_serovar": "app_serovar"}
        )[["sample", "app_serovar"]]

        merged = suis_df.merge(app_clean, on="sample", how="outer")
        merged_out = results_dir / "combined_summary.tsv"
        merged.to_csv(merged_out, sep="\t", index=False)
        log(f"[SUCCESS] Combined summary written → {merged_out}")

if __name__ == "__main__":
    main()
