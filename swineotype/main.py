#!/usr/bin/env python3
"""
swineotype.py — Serotyping tool for S. suis (and now A. pleuropneumoniae)
"""

from __future__ import annotations

import os
import sys
import glob
import click
from pathlib import Path

from swineotype.stages import stage1_score, stage2_resolver_call, interpret_resolver
from swineotype.config import load_config
from swineotype.adapters.app import run_app_analysis
from swineotype.utils import ensure_tool, ensure_unix_line_endings

# -------- Main orchestration --------

def process_one(assembly: str, out_dir: Path, threads: int, config: dict):
    run_dir = out_dir / Path(assembly).stem; run_dir.mkdir(parents=True, exist_ok=True)
    assembly = ensure_unix_line_endings(assembly, config["tmp_dir"])
    s1 = stage1_score(assembly, config["wzxwzy_fasta"], threads, run_dir, config)
    s1_top, s1_second = s1.get("top"), s1.get("second")
    allowed_pair = None
    if (s1_top in config["pair_1_14"]) or (s1_second in config["pair_1_14"]): allowed_pair = "1_vs_14"
    elif (s1_top in config["pair_2_1_2"]) or (s1_second in config["pair_2_1_2"]): allowed_pair = "2_vs_1_2"
    must_stage2 = (not s1.get("decisive", False)) or s1.get("must_stage2_for_pair", False)
    s2_ev, s2_status = None, "SKIPPED"
    if must_stage2 and allowed_pair:
        s2_ev = stage2_resolver_call(assembly, config["resolver_refs_fasta"], threads, run_dir, config, allowed_pair)
        s2_status = "OK" if s2_ev else "NO_HSP_OR_LOW_QUAL"
    final_sero, final_status = None, None
    if s2_ev:
        final_sero = interpret_resolver(s2_ev, config); final_status = "STAGE2" if final_sero else "NO_CALL_STAGE2"
    elif s1.get("decisive", False) and not must_stage2:
        final_sero = s1.get("top"); final_status = "STAGE1"
    else:
        final_status = "NO_CALL_STAGE2"
    return {"sample":assembly,"stage1_top":s1.get("top") or "","ref_id":(s2_ev or {}).get("ref_id",""),
            "contig":(s2_ev or {}).get("contig",""),"contig_pos":(s2_ev or {}).get("contig_pos",""),
            "strand":(s2_ev or {}).get("strand",""),"base":(s2_ev or {}).get("base",""),
            "status":final_status,"final_serotype":final_sero or ""}

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
    for tool in ("blastn","makeblastdb","samtools"): ensure_tool(tool)
    assemblies = expand_globs(list(assembly))
    merged_rows = []
    with click.progressbar(assemblies, label="Serotyping assemblies") as bar:
        for asm in bar:
            row = process_one(asm,out_dir,threads, config); merged_rows.append(row)
            fname, status, final = Path(asm).name,row["status"],row["final_serotype"]
            if status in ("STAGE1","STAGE2"): click.echo(f"[OK] {fname} => {final} ({status})")
            else: click.echo(f"[WARN] {fname} => {status}", err=True)
    if merged_csv:
        mpath = Path(merged_csv); mpath.parent.mkdir(parents=True, exist_ok=True)
        if not mpath.exists():
            mpath.write_text("sample,stage1_top,ref_id,contig,contig_pos,strand,base,status,final_serotype\n")
        with mpath.open("a") as fh:
            for r in merged_rows:
                fh.write(",".join(str(r.get(k,"")) for k in ["sample","stage1_top","ref_id","contig","contig_pos","strand","base","status","final_serotype"])+"\n")
        click.echo(f"[INFO] Merged CSV written: {mpath}")

if __name__=="__main__": main()
