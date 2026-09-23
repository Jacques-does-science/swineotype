#!/usr/bin/env python3
"""Check and regenerate the Stage-2 resolver references.

The four resolver references are slices of four public records. Two of the
shipped sequences carried a two-base deletion relative to their source record,
which put the downstream two thirds of each gene out of frame and -- because
the deletion sat between the start codon and the diagnostic site -- left the
declared diagnostic position pointing at the right base only by accident.

Everything this script needs in order to *verify* the shipped file is pinned in
``data/suis_resolver_refs.manifest.json``: accession.version, source
coordinates, strand, diagnostic position and the SHA-256 of the normalized
sequence. ``check`` is therefore fully offline; only ``regenerate`` touches
the network, and it is never run by the test suite.

Commands
--------
  check       verify data/suis_resolver_refs.fasta against the manifest
  regenerate  re-extract from ENA and rewrite the FASTA (network; opt-in)

Normalized sequence
-------------------
The checksummed string is the uppercase, concatenated nucleotide sequence
encoded as ASCII, with no header, no whitespace and no trailing newline.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
FASTA_PATH = REPO_ROOT / "data" / "suis_resolver_refs.fasta"
MANIFEST_PATH = REPO_ROOT / "data" / "suis_resolver_refs.manifest.json"

ENA_URL = "https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=false"

LINE_WIDTH = 60

STOP_CODONS = {"TAA", "TAG", "TGA"}

# --- normalized sequence helpers ---------------------------------------


def normalized(seq: str) -> str:
    """The exact string the manifest checksums: uppercase, no whitespace."""
    return "".join(seq.split()).upper()


def sha256_of(seq: str) -> str:
    return hashlib.sha256(normalized(seq).encode("ascii")).hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """[(header, sequence)] in file order. Headers keep their `|`-tag suffix."""
    records: list[tuple[str, list[str]]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                records.append((line[1:], []))
            elif line and records:
                records[-1][1].append(line)
    return [(h, normalized("".join(parts))) for h, parts in records]


def write_fasta(path: Path, records) -> Path:
    """Write a mapping, or (header, sequence) pairs, as wrapped FASTA."""
    with open(path, "w") as fh:
        for header, seq in dict(records).items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), LINE_WIDTH):
                fh.write(seq[i:i + LINE_WIDTH] + "\n")
    return path


def ref_id(header: str) -> str:
    return header.split("|", 1)[0]


# --- coding-sequence integrity -----------------------------------------


def coding_problems(seq: str) -> list[str]:
    """Ways in which ``seq`` fails to be a single intact coding sequence."""
    problems = []
    if len(seq) % 3:
        problems.append(f"length {len(seq)} is not a multiple of 3")
    if seq[:3] != "ATG":
        problems.append(f"start codon is {seq[:3]!r}, not ATG")
    if seq[-3:] not in STOP_CODONS:
        problems.append(f"final codon is {seq[-3:]!r}, not a stop codon")
    internal = [i // 3 + 1 for i in range(0, len(seq) - 3, 3)
                if seq[i:i + 3] in STOP_CODONS]
    if internal:
        problems.append(f"internal stop codon(s) at codon {internal}")
    return problems


# --- manifest ----------------------------------------------------------


def load_manifest(path: Path = MANIFEST_PATH) -> dict:
    with open(path) as fh:
        return json.load(fh)


def check(fasta_path: Path = FASTA_PATH, manifest_path: Path = MANIFEST_PATH) -> list[str]:
    """Verify the shipped FASTA against the pinned manifest. Offline.

    Returns a list of human-readable failures; empty means the file is exactly
    the pinned reference set.
    """
    manifest = load_manifest(manifest_path)
    entries = {e["id"]: e for e in manifest["references"]}
    records = {ref_id(h): (h, s) for h, s in read_fasta(fasta_path)}

    failures = []
    if set(records) != set(entries):
        failures.append(f"reference ids {sorted(records)} != manifest {sorted(entries)}")

    for name, entry in entries.items():
        if name not in records:
            continue
        header, seq = records[name]
        if header != entry["fasta_header"]:
            failures.append(f"{name}: header {header!r} != manifest {entry['fasta_header']!r}")
        if len(seq) != entry["length"]:
            failures.append(f"{name}: length {len(seq)} != manifest {entry['length']}")
        digest = sha256_of(seq)
        if digest != entry["sha256"]:
            failures.append(f"{name}: sha256 {digest} != manifest {entry['sha256']}")
        pos = entry["diagnostic_position"]
        if len(seq) >= pos:
            codon = seq[pos - 3:pos]
            if codon != entry["diagnostic_codon"]:
                failures.append(f"{name}: codon at {pos} is {codon!r} != manifest "
                                f"{entry['diagnostic_codon']!r}")
            if pos % 3:
                failures.append(f"{name}: diagnostic position {pos} is not a codon-third base")
        else:
            failures.append(f"{name}: sequence is shorter than diagnostic position {pos}")
        for problem in coding_problems(seq):
            failures.append(f"{name}: {problem}")
        expected_span = entry["source_end"] - entry["source_start"] + 1
        if expected_span != entry["length"]:
            failures.append(f"{name}: manifest coordinates span {expected_span} "
                            f"but declare length {entry['length']}")
    return failures


# --- regeneration from the public records ------------------------------


def fetch_embl(accession: str) -> str:
    """Download one EMBL flatfile from ENA. Network; never called by tests."""
    from urllib.request import urlopen
    url = ENA_URL.format(accession=accession)
    with urlopen(url, timeout=120) as resp:
        return resp.read().decode("ascii", errors="replace")


def sequence_from_embl(text: str) -> str:
    """The nucleotide sequence of an EMBL flatfile, uppercased and unspaced."""
    lines = text.splitlines()
    try:
        start = next(i for i, ln in enumerate(lines) if ln.startswith("SQ   "))
    except StopIteration:
        raise ValueError("no SQ line in EMBL record")
    chunks = []
    for ln in lines[start + 1:]:
        if ln.startswith("//"):
            break
        # Sequence lines are indented and end with a running coordinate.
        chunks.append("".join(ch for ch in ln if ch.isalpha()))
    return "".join(chunks).upper()


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def regenerate(fasta_path: Path = FASTA_PATH) -> list[str]:
    """Re-extract every reference from ENA and rewrite the FASTA.

    Verifies each extracted slice against the manifest checksum before writing,
    so a changed upstream record fails loudly instead of silently reshipping a
    different sequence.
    """
    manifest = load_manifest()
    notes, records = [], []
    for entry in manifest["references"]:
        embl = fetch_embl(entry["accession"])
        full = sequence_from_embl(embl)
        sub = full[entry["source_start"] - 1:entry["source_end"]]
        if entry["strand"] == "-":
            sub = reverse_complement(sub)
        digest = sha256_of(sub)
        if digest != entry["sha256"]:
            raise SystemExit(
                f"[ERROR] {entry['id']}: extracted sha256 {digest} != manifest "
                f"{entry['sha256']}. The upstream record may have changed; "
                f"investigate before updating the manifest.")
        problems = coding_problems(sub)
        if problems:
            raise SystemExit(f"[ERROR] {entry['id']}: extracted sequence is not an intact CDS: {problems}")
        notes.append(f"{entry['id']}: {entry['accession']}:"
                     f"{entry['source_start']}-{entry['source_end']}({entry['strand']}) OK")
        records.append((entry["fasta_header"], sub))
    write_fasta(fasta_path, records)
    return notes


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["check", "regenerate"])
    parser.add_argument("--fasta", type=Path, default=FASTA_PATH)
    args = parser.parse_args(argv)

    if args.command == "check":
        failures = check(args.fasta)
        for f in failures:
            print(f"[FAIL] {f}")
        if failures:
            return 1
        print(f"[OK] {args.fasta} matches {MANIFEST_PATH.name}")
        return 0

    notes = regenerate(args.fasta)
    for n in notes:
        print(f"[OK] {n}")
    failures = check(args.fasta)
    for f in failures:
        print(f"[FAIL] {f}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
