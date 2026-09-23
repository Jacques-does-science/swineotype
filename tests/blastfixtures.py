"""Builders for real-BLAST fixtures derived from the shipped reference data.

These assemble synthetic contigs out of the project's own references. They
exercise the real blastn/makeblastdb path, so they say whether the pipeline
behaves correctly on sequences of this shape -- they say nothing about how
often such sequences occur, and no accuracy figure may be read off them.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from resolver_refs import load_manifest, read_fasta, sha256_of, write_fasta  # via conftest.py

from swineotype.config import load_config

# Where the two-base (TG) deletion sat in the references as they shipped
# before 0.2.0, so the same fixture can be run against that damaged set.
DAMAGE_AT = {"cps2K": 881, "cps14K": 890}


def resolver_records() -> dict[str, str]:
    return dict(read_fasta(load_config()["resolver_refs_fasta"]))


def whitelist_records() -> dict[str, str]:
    return dict(read_fasta(load_config()["wzxwzy_fasta"]))


def by_id(records: dict[str, str]) -> dict[str, str]:
    return {h.split("|")[0].split()[0]: s for h, s in records.items()}


def damaged_resolver_fasta(path: Path) -> Path:
    """The reference set as it shipped, two bases short in cps2K and cps14K."""
    shipped = {e["id"]: e["damaged_sha256"] for e in load_manifest()["references"]}
    out = {}
    for header, seq in resolver_records().items():
        name = header.split("|")[0]
        if name in DAMAGE_AT:
            at = DAMAGE_AT[name]
            seq = seq[:at] + seq[at + 2:]
            assert sha256_of(seq) == shipped[name], f"{name}: not the state that shipped"
        out[header] = seq
    return write_fasta(path, out)


def markers_from(accession: str) -> dict[str, str]:
    """The shipped wzx/wzy references drawn from one source record."""
    found = {h.split()[0]: s for h, s in whitelist_records().items()
             if h.split()[0].endswith(accession)}
    assert found, f"no marker references from {accession}"
    return found


def type2_markers() -> dict[str, str]:
    """wzx and wzy of serotype 2, straight from the shipped marker panel."""
    return markers_from("BR001000")


def cps2k_sequence() -> str:
    return by_id(resolver_records())["cps2K"]


def with_base_at(seq: str, index0: int, base: str) -> str:
    """`seq` with sequence index `index0` (zero-based) replaced by `base`."""
    return seq[:index0] + base + seq[index0 + 1:]


def spacer(n: int = 500) -> str:
    """Filler so contigs are not flush against a gene boundary."""
    return "AGCT" * (n // 4)


def build_assembly(path: Path, resolver_copies: dict[str, str],
                   include_markers: bool = True,
                   markers: dict[str, str] | None = None) -> Path:
    """An assembly carrying wzx/wzy markers plus named cpsK copies.

    Each cpsK copy goes on its own contig, which is what makes them distinct
    physical loci; the same sequence repeated on one contig would be one locus.
    """
    contigs: dict[str, str] = {}
    if include_markers:
        markers = markers if markers is not None else type2_markers()
        contigs["markers"] = spacer() + spacer().join(markers.values()) + spacer()
    for name, seq in resolver_copies.items():
        contigs[name] = spacer() + seq + spacer()
    return write_fasta(path, contigs)


def blast_config(tmp_path: Path, **over) -> dict:
    cfg = load_config()
    cfg["tmp_dir"] = tmp_path / "cache"
    cfg["tmp_dir"].mkdir(parents=True, exist_ok=True)
    cfg["keep_debug"] = 0
    cfg["gzip_debug"] = 0
    cfg.update(over)
    return cfg


requires_blast = pytest.mark.skipif(
    not all(__import__("shutil").which(t) for t in ("blastn", "makeblastdb")),
    reason="blastn/makeblastdb not on PATH: real-BLAST integration not verified")
