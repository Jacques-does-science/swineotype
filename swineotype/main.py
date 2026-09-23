#!/usr/bin/env python3
"""
swineotype.py — Serotyping tool for S. suis (and now A. pleuropneumoniae)
"""

from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import json
import os
import sys
import glob
import click
from collections import Counter
from pathlib import Path

from swineotype import __version__
from swineotype.stages import (
    RESOLVABLE_FAMILIES,
    S2_SKIPPED,
    family_label,
    interpret_resolver,
    pair_of,
    resolver_status,
    stage1_score,
    stage2_resolver_call,
)
from swineotype.config import load_config
from swineotype.adapters.app import run_app_analysis
from swineotype.utils import ensure_tool, ensure_unix_line_endings, unique_run_names

SUMMARY_COLUMNS = ["sample", "sample_path", "run_dir",
                   "matched_reference_taxon", "input_species", "species_assessment",
                   "species",
                   "stage1_top", "stage1_family", "stage1_family_label",
                   "stage1_wzx_only", "stage2_status", "ref_id", "contig",
                   "contig_pos", "strand", "base", "triplet", "triplet_status",
                   "coding_status", "resolver_loci",
                   "status", "family_serotype", "final_serotype", "warnings"]

DEFAULT_SUMMARY_CSV = "swineotype_summary.csv"
DEFAULT_RUN_JSON = "swineotype_run.json"

# --- statuses ---
NO_CPS_MATCH = "NO_CPS_MATCH"
NO_WZY_MATCH = "NO_WZY_MATCH"
NON_TARGET_SPECIES = "NON_TARGET_SPECIES"
STAGE1 = "STAGE1"
STAGE2 = "STAGE2"
FAMILY_ONLY = "FAMILY_ONLY"
NO_CALL_FAMILY_AMBIGUOUS = "NO_CALL_FAMILY_AMBIGUOUS"
NO_CALL_PAIR_CONFLICT = "NO_CALL_PAIR_CONFLICT"

CALLED_STATUSES = (STAGE1, STAGE2)

# --- species assessment ---
SPECIES_NOT_ASSESSED = "NOT_ASSESSED"
SPECIES_USER_SUPPLIED = "USER_SUPPLIED"

# -------- Main orchestration --------

def family_view(s1: dict, config: dict) -> dict:
    """Family-level reading of a Stage-1 result.

    Confidence is established between FAMILIES -- {1,14}, {2,1/2}, and every
    other type on its own -- before any within-family resolution is attempted.
    Ambiguity inside a family is expected and is what Stage 2 settles;
    ambiguity between families is not, and cannot be settled by a site that
    only separates members of one family.

    A Stage-1 result carrying no family fields has not been assessed at family
    level, so it is reported as not decisive rather than assumed confident.
    """
    top = s1.get("family_top")
    return {"top": top,
            "second": s1.get("family_second"),
            "decisive": bool(s1.get("family_decisive", False)),
            "fraction": s1.get("family_fraction", 0.0),
            "delta": s1.get("family_delta", 0.0),
            "label": family_label(top, config)}


def _render_loci(s2_ev: dict | None) -> str:
    """Every distinct physical resolver locus, as one auditable field."""
    if not s2_ev:
        return ""
    out = []
    for ev in s2_ev.get("loci", [s2_ev]):
        out.append(f"{ev.get('locus', ev.get('contig', ''))}="
                   f"{ev.get('triplet', '')}"
                   f"[{ev.get('triplet_status', '')}]"
                   f"->{ev.get('implied_serotype') or '?'}")
    return ";".join(out)


def result_row(source, run_dir, s1, s2_ev, final_sero, final_status,
               matched_taxon, warnings, stage2_status=S2_SKIPPED,
               family_sero="", input_species=None, species_assessment=SPECIES_NOT_ASSESSED,
               fam=None):
    """One summary row.

    `sample` is the bare stem so it matches the key the APP adapter writes;
    `sample_path` keeps the provenance.

    Three separate columns carry what used to be conflated into `species`:
    `matched_reference_taxon` is the organism the matched *reference* is
    labelled with, `input_species` is an independently established identity if
    one was supplied, and `species_assessment` says which of those happened.
    `species` is retained as a deprecated alias of the first.
    """
    s2_ev = s2_ev or {}
    fam = fam or {}
    return {"sample": Path(source).stem, "sample_path": str(Path(source).resolve()),
            "run_dir": run_dir.name,
            "matched_reference_taxon": matched_taxon or "",
            "input_species": input_species or "",
            "species_assessment": species_assessment,
            "species": matched_taxon or "",
            "stage2_status": stage2_status,
            "stage1_top": (s1 or {}).get("top") or "",
            "stage1_family": fam.get("top") or "",
            "stage1_family_label": fam.get("label") or "",
            "stage1_wzx_only": (s1 or {}).get("top_wzx_only") or "",
            "ref_id": s2_ev.get("ref_id", ""),
            "contig": s2_ev.get("contig", ""), "contig_pos": s2_ev.get("contig_pos", ""),
            "strand": s2_ev.get("strand", ""), "base": s2_ev.get("base", ""),
            "triplet": s2_ev.get("triplet", ""),
            "triplet_status": s2_ev.get("triplet_status", ""),
            "coding_status": s2_ev.get("coding_status", ""),
            "resolver_loci": _render_loci(s2_ev or None),
            "status": final_status,
            "family_serotype": family_sero or "",
            "final_serotype": final_sero or "",
            "warnings": ";".join(warnings)}


def process_one(assembly: str, out_dir: Path, threads: int, config: dict,
                run_name: str | None = None, input_species: str | None = None):
    source = assembly
    run_dir = out_dir / (run_name or Path(assembly).stem); run_dir.mkdir(parents=True, exist_ok=True)
    assembly = ensure_unix_line_endings(assembly, config["tmp_dir"])
    warnings: list[str] = []
    s1 = stage1_score(assembly, config["wzxwzy_fasta"], threads, run_dir, config)
    s1_top = s1.get("top")
    fam = family_view(s1, config)

    species_assessment = SPECIES_USER_SUPPLIED if input_species else SPECIES_NOT_ASSESSED

    def row(status, *, final_sero=None, s2_ev=None, matched_taxon=None,
            warnings=(), stage2_status=S2_SKIPPED, family_sero=""):
        return result_row(source, run_dir, s1, s2_ev, final_sero, status,
                          matched_taxon, list(warnings),
                          stage2_status=stage2_status, family_sero=family_sero,
                          input_species=input_species,
                          species_assessment=species_assessment, fam=fam)

    # Nothing in the panel matched. This used to report the target species
    # anyway, because `species or config["target_species"]` filled the column
    # from a default: an assembly of an unrelated organism, or an empty file,
    # came back as "Streptococcus suis". Reference metadata is not an
    # identification of the input, and no match is not an identification at all.
    if s1_top is None and not s1.get("top_wzx_only"):
        return row(NO_CPS_MATCH, matched_taxon=None,
                   warnings=["no_cps_reference_matched"])

    # cps genes are present but no reference *wzy* matched. That establishes
    # only that no qualifying reference wzy was found -- NOT that the isolate
    # carries an intact locus, and not that the locus is novel. Report the
    # nearest wzx relative as a lead and say plainly what was not assessed.
    if s1_top is None:
        nearest = s1["top_wzx_only"]
        return row(NO_WZY_MATCH, matched_taxon=None,
                   warnings=["no_reference_wzy_matched",
                             f"nearest_wzx_relative:cps_type_{nearest}",
                             "locus_integrity_not_assessed"])

    # Species gate. The reference panel still carries the cps loci of six
    # former S. suis serotypes that have since been moved to other taxa, so a
    # best hit to one of those identifies a different organism. Reporting it as
    # an S. suis serotype -- which is what happened before -- is wrong at the
    # species level, not just the type level.
    matched_taxon = s1.get("top_species")
    if matched_taxon and matched_taxon != config["target_species"]:
        return row(NON_TARGET_SPECIES, matched_taxon=matched_taxon,
                   warnings=[f"cps_type_{s1_top}_belongs_to_{matched_taxon.replace(' ', '_')}"])

    # Family-level confidence, before any within-family resolution.
    if not fam["decisive"]:
        detail = f"{fam['top']}/{fam['second']}" if fam["second"] else f"{fam['top']}"
        warnings.append(f"competing_cps_families:{detail}")
        warnings.append(f"family_fraction={fam['fraction']:.2f};family_delta={fam['delta']:.0f}")
        return row(NO_CALL_FAMILY_AMBIGUOUS, matched_taxon=matched_taxon,
                   warnings=warnings)

    # A singleton family: the type IS the answer, no within-family site to read.
    # The runner-up is never a fallback into a resolvable family: with top="9"
    # and second="2" the tool used to read the cpsK site and report a
    # confident "2", from a site that only separates 2 from 1/2.
    if fam["top"] not in RESOLVABLE_FAMILIES:
        return row(STAGE1, final_sero=fam["top"].removeprefix("type:"),
                   matched_taxon=matched_taxon, warnings=warnings)

    # The FAMILY decision picks the resolver pair, not the top individual
    # label. The two can legitimately differ -- a family whose members each
    # contribute a different marker class can out-score the single highest
    # individual type -- and when they do, the family-level evidence is what
    # was actually assessed for confidence. Disagreement is worth recording.
    allowed_pair = fam["top"]
    top_label_pair = pair_of(s1_top, config)
    if top_label_pair and top_label_pair != allowed_pair:
        warnings.append(f"family_disagrees_with_top_label:{allowed_pair}/{top_label_pair}")

    s2_ev = stage2_resolver_call(assembly, config["resolver_refs_fasta"], threads,
                                 run_dir, config, allowed_pair)
    s2_status = resolver_status(s2_ev)

    if s2_ev:
        if s2_ev.get("conflict"):
            warnings.append(f"conflicting_resolver_copies:{','.join(s2_ev['conflicting_readings'])}")
            warnings.append(f"resolver_loci:{s2_ev.get('n_loci')}")
        elif s2_ev.get("triplet_status") != "OK":
            warnings.append(f"resolver_triplet_{s2_ev['triplet_status'].lower()}:{s2_ev.get('triplet')}")
        if s2_ev.get("coding_status") == "DISRUPTED":
            warnings.append(f"resolver_coding_disrupted:{s2_ev.get('coding_detail')}")
        elif s2_ev.get("coding_status") == "UNASSESSED":
            warnings.append("resolver_coding_integrity_unassessed")

    final_sero = interpret_resolver(s2_ev, config)

    # A Stage-2 call must stay inside the family Stage 1 pointed at.
    if final_sero and allowed_pair and pair_of(final_sero, config) not in (None, allowed_pair):
        warnings.append(f"stage2_outside_stage1_pair:{final_sero}")
        return row(NO_CALL_PAIR_CONFLICT, s2_ev=s2_ev, matched_taxon=matched_taxon,
                   warnings=warnings, stage2_status=s2_status, family_sero=fam["label"])

    if final_sero:
        return row(STAGE2, final_sero=final_sero, s2_ev=s2_ev,
                   matched_taxon=matched_taxon, warnings=warnings, stage2_status=s2_status)

    # The family is established even though the exact member is not. Report it
    # rather than throwing the whole result away.
    return row(FAMILY_ONLY, s2_ev=s2_ev, matched_taxon=matched_taxon,
               warnings=warnings, stage2_status=s2_status, family_sero=fam["label"])

# -------- CLI --------

def expand_globs(paths: list[str]) -> tuple[list[str], list[str]]:
    """Expand glob patterns. Returns (matched paths, patterns that matched nothing).

    An unmatched pattern used to expand to nothing and the run then finished
    successfully with an empty summary -- a typo'd path was indistinguishable
    from a batch where every sample was uninformative.
    """
    out: list[str] = []
    unmatched: list[str] = []
    for p in paths:
        if any(ch in p for ch in "*?["):
            matches = sorted(glob.glob(p))
            if not matches:
                unmatched.append(p)
            out.extend(matches)
        elif not Path(p).exists():
            unmatched.append(p)
        else:
            out.append(p)
    return out, unmatched


def _file_digest(path) -> str | None:
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return None


def run_record(config: dict, assemblies: list[str], rows: list[dict],
               input_species: str | None) -> dict:
    """Everything needed to reproduce and audit this run."""
    manifest_version = None
    manifest_path = Path(config["data_dir"]) / "suis_resolver_refs.manifest.json"
    try:
        manifest_version = json.loads(manifest_path.read_text()).get("data_version")
    except (OSError, ValueError):
        pass
    return {
        "swineotype_version": __version__,
        "run_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "argv": sys.argv,
        "python": sys.version.split()[0],
        "input_species": input_species or "",
        "assemblies": [str(Path(a).resolve()) for a in assemblies],
        "reference_data": {
            "resolver_refs_fasta": str(config["resolver_refs_fasta"]),
            "resolver_refs_sha256": _file_digest(config["resolver_refs_fasta"]),
            "resolver_refs_manifest_version": manifest_version,
            "wzxwzy_fasta": str(config["wzxwzy_fasta"]),
            "wzxwzy_sha256": _file_digest(config["wzxwzy_fasta"]),
        },
        "effective_config": config,
        "status_counts": dict(Counter(r["status"] for r in rows)),
    }


def write_summary_csv(path: Path, rows: list[dict], append: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not (append and path.exists())

    if append and path.exists():
        # Appending new-shaped rows under an old header writes every value
        # into the wrong column, silently. The summary gained and lost columns
        # in 0.2.0, so a merged CSV from 0.1.x must not be appended to.
        with path.open(newline="") as fh:
            existing = next(csv.reader(fh), None)
        if existing is not None and existing != SUMMARY_COLUMNS:
            raise click.ClickException(
                f"{path} has a different column set (probably written by an "
                f"older swineotype). Appending would misalign every value. "
                f"Use a new --merged_csv path, or move the old file aside.\n"
                f"  existing: {existing}\n  expected: {SUMMARY_COLUMNS}")

    mode = "a" if append else "w"
    # csv.writer so a comma or quote in a path cannot corrupt the row.
    with path.open(mode, newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        for r in rows:
            writer.writerow(r)


@click.command()
@click.option("--assembly", multiple=True, required=True, type=click.Path(), help="Path to one or more assembly files. Globs are supported.")
@click.option("--out_dir", required=True, type=click.Path(), help="Output directory")
@click.option("--merged_csv", default=None, type=click.Path(), help="Additional CSV to append all results to. A per-run summary is always written to the output directory regardless.")
@click.option("--threads", default=lambda: max(1, os.cpu_count() // 2), help="Number of threads to use")
@click.option("--species", default="suis", type=click.Choice(["suis", "app"]), help="Species to serotype")
@click.option("--input_species", default=None, help="Independently established species of the input (e.g. from ANI). Recorded as-is; this tool does not measure it.")
@click.option("--config", default=None, type=click.Path(exists=True), help="Path to a custom config.yaml file")
@click.version_option(__version__, prog_name="swineotype")
def main(assembly, out_dir, merged_csv, threads, species, input_species, config):
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

    assemblies, unmatched = expand_globs(list(assembly))
    if unmatched:
        for p in unmatched:
            click.echo(f"[ERROR] No assembly matched: {p}", err=True)
        raise SystemExit(2)
    if not assemblies:
        click.echo("[ERROR] No assemblies to process", err=True)
        raise SystemExit(2)

    out_dir = Path(out_dir).resolve(); out_dir.mkdir(parents=True, exist_ok=True)
    # Keep staged assemblies and BLAST databases with the run's own outputs.
    # They used to land in <install>/data/tmp, which persists across unrelated
    # runs and breaks outright on a read-only (e.g. shared conda) install.
    if not config.get("tmp_dir_explicit"):
        config["tmp_dir"] = out_dir / ".swineotype_cache"
        config["tmp_dir"].mkdir(parents=True, exist_ok=True)
    # samtools is no longer used: Stage 2 reads the diagnostic codon out of the
    # BLAST alignment itself.
    for tool in ("blastn","makeblastdb"): ensure_tool(tool)
    run_names = unique_run_names(assemblies)
    merged_rows = []
    with click.progressbar(list(zip(assemblies, run_names)), label="Serotyping assemblies") as bar:
        for asm, run_name in bar:
            row = process_one(asm, out_dir, threads, config, run_name,
                              input_species=input_species)
            merged_rows.append(row)
            fname, status, final = Path(asm).name,row["status"],row["final_serotype"]
            if status in CALLED_STATUSES: click.echo(f"[OK] {fname} => {final} ({status})")
            elif status == FAMILY_ONLY:
                click.echo(f"[WARN] {fname} => {row['family_serotype']} "
                           f"(family only; exact resolution withheld: {row['stage2_status']})", err=True)
            elif status == NO_CPS_MATCH:
                click.echo(f"[WARN] {fname} => no cps reference matched; no species inferred", err=True)
            elif status == NO_WZY_MATCH:
                click.echo(f"[WARN] {fname} => no reference wzy matched "
                           f"(nearest wzx: cps type {row['stage1_wzx_only']}); "
                           f"locus integrity not assessed", err=True)
            elif status == NON_TARGET_SPECIES:
                click.echo(f"[WARN] {fname} => matched reference is not {config['target_species']}: "
                           f"{row['matched_reference_taxon']} (cps type {row['stage1_top']})", err=True)
            else: click.echo(f"[WARN] {fname} => {status}", err=True)

    # Structured results and the effective configuration are persisted by
    # default. They used to exist only if --merged_csv was passed, so the
    # common invocation left nothing machine-readable and nothing recording
    # which thresholds and reference data produced the calls.
    summary_path = out_dir / DEFAULT_SUMMARY_CSV
    write_summary_csv(summary_path, merged_rows, append=False)
    record = run_record(config, assemblies, merged_rows, input_species)
    # json.dumps recurses on its own; only the leaves it cannot encode need
    # help. Sets are sorted so the record is stable between runs.
    (out_dir / DEFAULT_RUN_JSON).write_text(json.dumps(
        record, indent=2,
        default=lambda o: sorted(o) if isinstance(o, (set, frozenset)) else str(o)) + "\n")
    click.echo(f"[INFO] Summary written: {summary_path}")
    click.echo(f"[INFO] Run record written: {out_dir / DEFAULT_RUN_JSON}")

    if merged_csv:
        mpath = Path(merged_csv)
        write_summary_csv(mpath, merged_rows, append=True)
        click.echo(f"[INFO] Merged CSV written: {mpath}")

if __name__=="__main__": main()
