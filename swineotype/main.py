#!/usr/bin/env python3
"""
swineotype.py — Serotyping tool for S. suis (and now A. pleuropneumoniae)
"""

from __future__ import annotations

import csv
import hashlib
import os
import sys
import glob
import click
from collections import Counter
from pathlib import Path

from swineotype.stages import stage1_score, stage2_resolver_call, interpret_resolver
from swineotype.config import load_config
from swineotype.adapters.app import run_app_analysis
from swineotype.utils import ensure_tool, ensure_unix_line_endings

SUMMARY_COLUMNS = ["sample", "sample_path", "run_dir", "species", "stage1_top",
                   "ref_id", "contig", "contig_pos", "strand", "base", "status",
                   "final_serotype", "warnings"]

# -------- Main orchestration --------

def pair_of(serotype: str | None, config: dict) -> str | None:
    """Which resolver pair a serotype belongs to, or None."""
    if serotype in config["pair_1_14"]: return "1_vs_14"
    if serotype in config["pair_2_1_2"]: return "2_vs_1_2"
    return None


def choose_pair(s1_top: str | None, s1_second: str | None, config: dict) -> str | None:
    """Pick the resolver pair implied by Stage 1.

    The top hit is checked against BOTH pairs before the runner-up is
    considered. The previous form tested pair_1_14 against top *and* second
    before ever testing pair_2_1_2, so a runner-up from the 1/14 group
    overrode a top hit from the 2/(1/2) group -- and because the four resolver
    references are 98-99.5% identical, Stage 2 would then happily return a
    confident serotype from the wrong pair.
    """
    return pair_of(s1_top, config) or pair_of(s1_second, config)


def unique_run_names(paths: list[str]) -> list[str]:
    """Per-sample output directory names, disambiguated only where needed.

    Two inputs sharing a basename would otherwise write their debug TSVs into
    the same directory and overwrite each other -- the same collision that
    affected the BLAST cache, and the one that makes a bad call impossible to
    investigate afterwards.
    """
    stems = [Path(p).stem for p in paths]
    duplicated = {s for s, n in Counter(stems).items() if n > 1}
    names = []
    for path, stem in zip(paths, stems):
        if stem in duplicated:
            digest = hashlib.sha1(str(Path(path).resolve()).encode()).hexdigest()[:8]
            names.append(f"{stem}__{digest}")
        else:
            names.append(stem)
    return names


def result_row(source, run_dir, s1, s2_ev, final_sero, final_status, species, warnings):
    """One summary row.

    `sample` is the bare stem so it matches the key the APP adapter writes;
    `sample_path` keeps the provenance. It previously held the *staged tmp*
    path, which meant the suis/APP merge on "sample" could never match.
    """
    s2_ev = s2_ev or {}
    return {"sample": Path(source).stem, "sample_path": str(Path(source).resolve()),
            "run_dir": run_dir.name, "species": species or "",
            "stage1_top": (s1 or {}).get("top") or "", "ref_id": s2_ev.get("ref_id", ""),
            "contig": s2_ev.get("contig", ""), "contig_pos": s2_ev.get("contig_pos", ""),
            "strand": s2_ev.get("strand", ""), "base": s2_ev.get("base", ""),
            "status": final_status, "final_serotype": final_sero or "",
            "warnings": ";".join(warnings)}


def process_one(assembly: str, out_dir: Path, threads: int, config: dict, run_name: str | None = None):
    source = assembly
    run_dir = out_dir / (run_name or Path(assembly).stem); run_dir.mkdir(parents=True, exist_ok=True)
    assembly = ensure_unix_line_endings(assembly, config["tmp_dir"])
    warnings: list[str] = []
    s1 = stage1_score(assembly, config["wzxwzy_fasta"], threads, run_dir, config)
    s1_top, s1_second = s1.get("top"), s1.get("second")

    # Species gate. The reference panel still carries the cps loci of six
    # former S. suis serotypes that have since been moved to other taxa, so a
    # best hit to one of those identifies a different organism. Reporting it as
    # an S. suis serotype -- which is what happened before -- is wrong at the
    # species level, not just the type level.
    species = s1.get("top_species") if s1_top else None
    if s1_top and species and species != config["target_species"]:
        return result_row(source, run_dir, s1, None, None,
                          "NON_TARGET_SPECIES", species,
                          [f"cps_type_{s1_top}_belongs_to_{species.replace(' ', '_')}"])

    allowed_pair = choose_pair(s1_top, s1_second, config)

    # Flag a cross-pair runner-up only when the top hit is not comfortably
    # ahead. Serotypes 1 and 14 share a locus and always draw serotype 2 into
    # second place, so warning on the configuration alone fires on every 1/14
    # isolate; it is the thin margin, not the disagreement, that is the risk.
    top_pair, second_pair = pair_of(s1_top, config), pair_of(s1_second, config)
    if top_pair and second_pair and top_pair != second_pair \
            and s1.get("delta", 0.0) < config["delta"]:
        warnings.append(f"stage1_pair_ambiguous:{s1_top}/{s1_second}")

    must_stage2 = (not s1.get("decisive", False)) or s1.get("must_stage2_for_pair", False)
    s2_ev, s2_status = None, "SKIPPED"
    if must_stage2 and allowed_pair:
        s2_ev = stage2_resolver_call(assembly, config["resolver_refs_fasta"], threads, run_dir, config, allowed_pair)
        s2_status = "OK" if s2_ev else "NO_HSP_OR_LOW_QUAL"
    if s2_ev and s2_ev.get("base") == "-":
        warnings.append("resolver_site_deleted")
    elif s2_ev and s2_ev.get("base") not in ("G", "C", "T"):
        warnings.append(f"resolver_base_not_gct:{s2_ev.get('base')}")

    final_sero, final_status = None, None
    if s2_ev:
        final_sero = interpret_resolver(s2_ev, config); final_status = "STAGE2" if final_sero else "NO_CALL_STAGE2"
    elif s1.get("decisive", False) and not must_stage2:
        final_sero = s1.get("top"); final_status = "STAGE1"
    else:
        final_status = "NO_CALL_STAGE2"

    # A Stage-2 call must stay inside the pair Stage 1 pointed at. This should
    # be unreachable now that choose_pair() prefers the top hit, but the cost of
    # a wrong answer here is a serotype 2 <-> 14 confusion, so it is checked.
    if final_sero and allowed_pair and pair_of(final_sero, config) not in (None, allowed_pair):
        warnings.append(f"stage2_outside_stage1_pair:{final_sero}")
        final_sero, final_status = None, "NO_CALL_PAIR_CONFLICT"

    return result_row(source, run_dir, s1, s2_ev, final_sero, final_status,
                      species or config["target_species"], warnings)

# -------- CLI --------

def expand_globs(paths: list[str]) -> list[str]:
    out = []
    for p in paths:
        if any(ch in p for ch in "*?[]"): out.extend(sorted(glob.glob(p)))
        else: out.append(p)
    return out

@click.command()
@click.option("--assembly", multiple=True, required=True, type=click.Path(), help="Path to one or more assembly files. Globs are supported.")
@click.option("--out_dir", required=True, type=click.Path(), help="Output directory")
@click.option("--merged_csv", default=None, type=click.Path(), help="Path to merge results into a single CSV file")
@click.option("--threads", default=lambda: max(1, os.cpu_count() // 2), help="Number of threads to use")
@click.option("--species", default="suis", type=click.Choice(["suis", "app"]), help="Species to serotype")
@click.option("--config", default=None, type=click.Path(exists=True), help="Path to a custom config.yaml file")
def main(assembly, out_dir, merged_csv, threads, species, config):
    """Swineotype: serotyping from assemblies"""
    config = load_config(config)

    if species=="app":
        run_app_analysis(
            assembly=list(assembly),
            out_dir=out_dir,
            threads=threads,
            swineotype_summary=merged_csv,
        )
        sys.exit(0)

    out_dir = Path(out_dir).resolve(); out_dir.mkdir(parents=True, exist_ok=True)
    # Keep staged assemblies and BLAST databases with the run's own outputs.
    # They used to land in <install>/data/tmp, which persists across unrelated
    # runs and breaks outright on a read-only (e.g. shared conda) install.
    if not config.get("tmp_dir_explicit"):
        config["tmp_dir"] = out_dir / ".swineotype_cache"
        config["tmp_dir"].mkdir(parents=True, exist_ok=True)
    # samtools is no longer used: Stage 2 reads the diagnostic base out of the
    # BLAST alignment itself.
    for tool in ("blastn","makeblastdb"): ensure_tool(tool)
    assemblies = expand_globs(list(assembly))
    run_names = unique_run_names(assemblies)
    merged_rows = []
    with click.progressbar(list(zip(assemblies, run_names)), label="Serotyping assemblies") as bar:
        for asm, run_name in bar:
            row = process_one(asm,out_dir,threads, config, run_name); merged_rows.append(row)
            fname, status, final = Path(asm).name,row["status"],row["final_serotype"]
            if status in ("STAGE1","STAGE2"): click.echo(f"[OK] {fname} => {final} ({status})")
            elif status == "NON_TARGET_SPECIES":
                click.echo(f"[WARN] {fname} => not {config['target_species']}: {row['species']} "
                           f"(cps type {row['stage1_top']})", err=True)
            else: click.echo(f"[WARN] {fname} => {status}", err=True)
    if merged_csv:
        mpath = Path(merged_csv); mpath.parent.mkdir(parents=True, exist_ok=True)
        write_header = not mpath.exists()
        # csv.writer so a comma or quote in a path cannot corrupt the row.
        with mpath.open("a", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS, extrasaction="ignore")
            if write_header: writer.writeheader()
            for r in merged_rows: writer.writerow(r)
        click.echo(f"[INFO] Merged CSV written: {mpath}")

if __name__=="__main__": main()
