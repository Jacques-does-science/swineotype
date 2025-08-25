#!/usr/bin/env python3
"""
swineotype.py — Serotyping tool for S. suis (and now A. pleuropneumoniae)

Stage-1: BLAST wzx/wzy whitelist against each assembly → sum bitscores per serotype.
          Decisiveness = plurality >= PLURALITY AND (top - second) >= DELTA.
          (Optionally require Wzx/Wzy agreement.)
          Force Stage-2 for ambiguous families even if decisive.

Stage-2: BLAST resolver references (e.g., cpsK for S. suis) → find HSPs that span the diagnostic site,
          but ONLY within the Stage-1 suggested family.
          Diagnostic site coordinates and expected alleles come from FASTA header metadata.

Outputs (per-sample folder):
  - wzxwzy_vs_asm.tsv  (Stage-1 BLAST table)
  - resolver_vs_asm.tsv (Stage-2 BLAST table)
  - swineotype.json     (rich summary)
  - swineotype.tsv      (compact one-line TSV)

Merged CSV (if requested):
  sample,stage1_top,ref_id,contig,contig_pos,strand,base,status,final_serotype
"""

from __future__ import annotations

import os
import sys
import json
import glob
import shlex
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict

# -------- Configuration --------

# Defaults (overridable by env or CLI species flag)
WZXWZY_FASTA = os.environ.get("SWINEO_WZXWZY", None)
RESOLVER_REFS_FASTA = os.environ.get("SWINEO_RESOLVER_REFS", None)

TMPDIR = Path(os.environ.get("SWINEO_TMPDIR", "results/db_cache")).resolve()
TMPDIR.mkdir(parents=True, exist_ok=True)

PLURALITY = float(os.environ.get("SWINEO_PLURALITY", "0.60"))
DELTA = float(os.environ.get("SWINEO_DELTA", "100"))
REQUIRE_AGREEMENT = int(os.environ.get("SWINEO_REQUIRE_AGREEMENT", "1"))

MIN_PID = float(os.environ.get("SWINEO_MIN_PID", "85.0"))
MIN_COV = float(os.environ.get("SWINEO_MIN_COV", "0.80"))

MIN_RES_PID  = float(os.environ.get("SWINEO_MIN_RESOLVER_PID", "90.0"))
MIN_RES_ALEN = int(os.environ.get("SWINEO_MIN_RESOLVER_ALEN", "300"))

KEEP_DEBUG   = int(os.environ.get("KEEP_DEBUG", "1"))
GZIP_DEBUG   = int(os.environ.get("GZIP_DEBUG", "0"))
CLEAN_TEMP   = int(os.environ.get("CLEAN_TEMP", "0"))

AMBIG_SET = {"1", "14", "2", "1/2"}
PAIR_1_14 = {"1", "14"}
PAIR_2_1_2 = {"2", "1/2"}

# -------- Utilities --------

def parse_resolver_meta(qid: str):
    """
    Parse resolver FASTA IDs of the form:
      >id|pair=1_vs_14|pos=483|baseA=G|A=1|B=14
    Returns dict with keys: pair, pos, baseA, A, B
    """
    meta = {"pair": None, "pos": None, "baseA": "G", "A": None, "B": None}
    for tok in qid.split("|")[1:]:
        if tok.startswith("pair="):   meta["pair"]   = tok.split("=",1)[1]
        elif tok.startswith("pos="):  meta["pos"]    = int(tok.split("=",1)[1])
        elif tok.startswith("baseA="):meta["baseA"]  = tok.split("=",1)[1].upper()
        elif tok.startswith("A="):    meta["A"]      = tok.split("=",1)[1]
        elif tok.startswith("B="):    meta["B"]      = tok.split("=",1)[1]
    return meta

def ensure_tool(name: str):
    from shutil import which
    if which(name) is None:
        sys.exit(f"[ERROR] Required tool not found in PATH: {name}")

def run(cmd, check=True, capture=True, cwd=None, text=True):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    res = subprocess.run(cmd, check=check, capture_output=capture, cwd=cwd, text=text)
    return res.stdout

def make_db_if_needed(asm_fa: str, tmpdir: Path) -> str:
    prefix = tmpdir / ("asmdb_" + Path(asm_fa).stem)
    nin = prefix.with_suffix(".nin")
    ndb = prefix.with_suffix(".ndb")
    if not (nin.exists() or ndb.exists()):
        run(["makeblastdb", "-in", asm_fa, "-dbtype", "nucl", "-out", str(prefix)])
    return str(prefix)

def run_blast(query_fa: str, db_prefix: str, threads: int, outfmt_cols: str, max_target_seqs=50) -> str:
    cmd = [
        "blastn",
        "-query", query_fa,
        "-db", db_prefix,
        "-task", "blastn",
        "-outfmt", outfmt_cols,
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
    ]
    return run(cmd)

def parse_whitelist_headers(fasta_path: str):
    allele_to_type = {}
    allele_to_geneclass = {}
    with open(fasta_path, "r") as fh:
        for line in fh:
            if not line.startswith(">"): continue
            h = line[1:].strip()
            allele_id = h.split()[0]
            st = None
            for tok in h.split():
                if tok.startswith("[type_id=") and tok.endswith("]"):
                    st = tok[len("[type_id="):-1]; break
            allele_to_type[allele_id] = st
            low = allele_id.lower()
            geneclass = "wzy" if "wzy" in low else ("wzx" if "wzx" in low else None)
            allele_to_geneclass[allele_id] = geneclass
    return allele_to_type, allele_to_geneclass

def infer_pair_from_ref(qid: str) -> str | None:
    s = qid.lower()
    if "1_vs_14" in s: return "1_vs_14"
    if "2_vs_1_2" in s: return "2_vs_1_2"
    return None

# -------- Stage-1 --------

def stage1_score(assembly_fa: str, whitelist_fa: str, threads: int, run_dir: Path):
    ensure_tool("blastn"); ensure_tool("makeblastdb")
    allele_to_type, allele_to_geneclass = parse_whitelist_headers(whitelist_fa)
    db_prefix = make_db_if_needed(assembly_fa, TMPDIR)
    outfmt = "6 qseqid sseqid pident length qlen evalue bitscore qstart qend sstart send"
    tsv_text = run_blast(whitelist_fa, db_prefix, threads, outfmt)
    stage1_tsv = run_dir / "wzxwzy_vs_asm.tsv"
    if KEEP_DEBUG:
        stage1_tsv.write_text(tsv_text + ("\n" if tsv_text else ""))
    best_by_q = {}
    score_by_type = defaultdict(float)
    for line in filter(None, tsv_text.splitlines()):
        qseqid, sseqid, pident, length, qlen, *_rest = line.split("\t")
        pident, length, qlen = float(pident), int(length), int(qlen)
        coverage = (length / qlen) if qlen else 0
        if pident < MIN_PID or coverage < MIN_COV: continue
        bitscore = float(_rest[3])
        st = allele_to_type.get(qseqid)
        if st: score_by_type[st] += bitscore
        if qseqid not in best_by_q or bitscore > best_by_q[qseqid]["bitscore"]:
            best_by_q[qseqid] = {"bitscore": bitscore}
    ordered = sorted(score_by_type.items(), key=lambda kv: kv[1], reverse=True)
    top, top_score = (ordered[0][0], ordered[0][1]) if ordered else (None, 0.0)
    second, second_score = (ordered[1][0], ordered[1][1]) if len(ordered) > 1 else (None, 0.0)
    total = sum(score_by_type.values())
    fraction, delta = (top_score/total if total else 0.0), top_score-second_score
    decisive = (fraction >= PLURALITY) and (delta >= DELTA)
    must_stage2_for_pair = bool(top in AMBIG_SET)
    return {"scores":score_by_type,"top":top,"second":second,"fraction":fraction,
            "delta":delta,"decisive":decisive,"must_stage2_for_pair":must_stage2_for_pair}

# -------- Stage-2 --------

def stage2_resolver_call(assembly_fa: str, resolver_refs_fa: str, threads: int, run_dir: Path, allowed_pair: str|None=None):
    ensure_tool("blastn"); ensure_tool("makeblastdb"); ensure_tool("samtools")
    db_prefix = make_db_if_needed(assembly_fa, TMPDIR)
    outfmt = "6 qseqid sseqid pident length qlen evalue bitscore qstart qend sstart send"
    tsv_text = run_blast(resolver_refs_fa, db_prefix, threads, outfmt)
    (run_dir / "resolver_vs_asm.tsv").write_text(tsv_text + ("\n" if tsv_text else ""))
    best = None
    for line in filter(None, tsv_text.splitlines()):
        qseqid, sseqid, pident, length, qlen, evalue, bitscore, qstart, qend, sstart, send = line.split("\t")
        pident, length, qstart, qend, sstart, send, bitscore = float(pident), int(length), int(qstart), int(qend), int(sstart), int(send), float(bitscore)
        meta = parse_resolver_meta(qseqid)
        if allowed_pair and meta["pair"] != allowed_pair: continue
        pos = meta["pos"]
        spans = (qstart <= pos <= qend) or (qend <= pos <= qstart)
        if not spans: continue
        if pident < MIN_RES_PID or length < MIN_RES_ALEN: continue
        qoff = pos - qstart
        if sstart <= send:
            strand = "+"; tpos = sstart + qoff
        else:
            strand = "-"; tpos = sstart - qoff
        ev = {"ref_id":qseqid,"contig":sseqid,"contig_pos":tpos,"strand":strand,
              "pident":pident,"length":length,"bitscore":bitscore,"pair":meta["pair"],"base":None}
        if best is None or bitscore > best[0]: best = (bitscore, ev)
    if not best: return None
    ev = best[1]
    region = f"{ev['contig']}:{ev['contig_pos']}-{ev['contig_pos']}"
    fa = run(["samtools","faidx",assembly_fa,region])
    lines = [ln.strip() for ln in fa.splitlines()]
    ev["base"] = lines[1].strip().upper() if len(lines)>1 else "N"
    return ev

def interpret_resolver(ev: dict|None) -> str|None:
    if ev is None: return None
    meta = parse_resolver_meta(ev["ref_id"])
    if meta["pair"] == "1_vs_14":
        return meta["B"] if ev["base"] == meta["baseA"] else meta["A"]
    if meta["pair"] == "2_vs_1_2":
        return meta["B"] if ev["base"] == meta["baseA"] else meta["A"]
    return None

# -------- Main orchestration --------

def process_one(assembly: str, out_dir: Path, threads: int):
    run_dir = out_dir / Path(assembly).stem; run_dir.mkdir(parents=True, exist_ok=True)
    s1 = stage1_score(assembly, WZXWZY_FASTA, threads, run_dir)
    s1_top, s1_second = s1.get("top"), s1.get("second")
    allowed_pair = None
    if (s1_top in PAIR_1_14) or (s1_second in PAIR_1_14): allowed_pair = "1_vs_14"
    elif (s1_top in PAIR_2_1_2) or (s1_second in PAIR_2_1_2): allowed_pair = "2_vs_1_2"
    must_stage2 = (not s1.get("decisive", False)) or s1.get("must_stage2_for_pair", False)
    s2_ev, s2_status = None, "SKIPPED"
    if must_stage2 and allowed_pair:
        s2_ev = stage2_resolver_call(assembly, RESOLVER_REFS_FASTA, threads, run_dir, allowed_pair)
        s2_status = "OK" if s2_ev else "NO_HSP_OR_LOW_QUAL"
    final_sero, final_status = None, None
    if s2_ev:
        final_sero = interpret_resolver(s2_ev); final_status = "STAGE2" if final_sero else "NO_CALL_STAGE2"
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

def main():
    parser = argparse.ArgumentParser(description="Swineotype: serotyping from assemblies")
    parser.add_argument("--assembly", nargs="+", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--merged_csv", required=False, default=None)
    parser.add_argument("--threads", type=int, default=max(1, os.cpu_count()//2))
    parser.add_argument("--species", choices=["suis","app"], default="suis")
    args = parser.parse_args()
    global WZXWZY_FASTA, RESOLVER_REFS_FASTA
    if args.species=="suis":
        WZXWZY_FASTA = WZXWZY_FASTA or "data/suis_wzxwzy_whitelist.fasta"
        RESOLVER_REFS_FASTA = RESOLVER_REFS_FASTA or "data/suis_resolver_refs.fasta"
    else:
        WZXWZY_FASTA = WZXWZY_FASTA or "data/app_wzxwzy_whitelist.fasta"
        RESOLVER_REFS_FASTA = RESOLVER_REFS_FASTA or "data/app_resolver_refs.fasta"
    out_dir = Path(args.out_dir).resolve(); out_dir.mkdir(parents=True, exist_ok=True)
    for tool in ("blastn","makeblastdb","samtools"): ensure_tool(tool)
    assemblies = expand_globs(args.assembly)
    merged_rows = []
    for asm in assemblies:
        row = process_one(asm,out_dir,args.threads); merged_rows.append(row)
        fname, status, final = Path(asm).name,row["status"],row["final_serotype"]
        if status in ("STAGE1","STAGE2"): print(f"[OK] {fname} => {final} ({status})")
        else: print(f"[WARN] {fname} => {status}")
    if args.merged_csv:
        mpath = Path(args.merged_csv); mpath.parent.mkdir(parents=True, exist_ok=True)
        if not mpath.exists():
            mpath.write_text("sample,stage1_top,ref_id,contig,contig_pos,strand,base,status,final_serotype\n")
        with mpath.open("a") as fh:
            for r in merged_rows:
                fh.write(",".join(str(r.get(k,"")) for k in ["sample","stage1_top","ref_id","contig","contig_pos","strand","base","status","final_serotype"])+"\n")
        print(f"[INFO] Merged CSV written: {mpath}")

if __name__=="__main__": main()
