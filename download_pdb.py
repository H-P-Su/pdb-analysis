#!/usr/bin/env python3
"""
Download PDB, CIF, and/or sequence files from RCSB (rcsb.org).

Usage:
    python download_pdb.py [OPTIONS] PDB_ID [PDB_ID ...]

Examples:
    python download_pdb.py 1ABC
    python download_pdb.py --cif 1ABC 2XYZ 3DEF
    python download_pdb.py --pdb --cif --seq 1ABC 2XYZ
    python download_pdb.py --seq --outdir ./sequences 1ABC
    python download_pdb.py --from-file ids.txt --pdb --cif
"""

import argparse
import sys
import urllib.request
import urllib.error
from pathlib import Path


RCSB_BASE = "https://files.rcsb.org/download"
RCSB_FASTA = "https://www.rcsb.org/fasta/entry"


def download_file(url: str, dest: Path, overwrite: bool = False) -> bool:
    if dest.exists() and not overwrite:
        print(f"  [skip] {dest.name} already exists (use --overwrite to re-download)")
        return True
    try:
        print(f"  Downloading {url} -> {dest.name}")
        urllib.request.urlretrieve(url, dest)
        return True
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"  [error] {dest.name}: not found (404) — check the PDB ID")
        else:
            print(f"  [error] {dest.name}: HTTP {e.code}")
        return False
    except urllib.error.URLError as e:
        print(f"  [error] {dest.name}: {e.reason}")
        return False


def download_entry(
    pdb_id: str,
    outdir: Path,
    fmt_pdb: bool,
    fmt_cif: bool,
    fmt_seq: bool,
    overwrite: bool,
) -> dict:
    pdb_id = pdb_id.upper().strip()
    results = {"id": pdb_id, "pdb": None, "cif": None, "seq": None}

    print(f"\n{pdb_id}")

    if fmt_pdb:
        url = f"{RCSB_BASE}/{pdb_id}.pdb"
        dest = outdir / f"{pdb_id}.pdb"
        results["pdb"] = download_file(url, dest, overwrite)

    if fmt_cif:
        url = f"{RCSB_BASE}/{pdb_id}.cif"
        dest = outdir / f"{pdb_id}.cif"
        results["cif"] = download_file(url, dest, overwrite)

    if fmt_seq:
        url = f"{RCSB_FASTA}/{pdb_id}"
        dest = outdir / f"{pdb_id}.fasta"
        results["seq"] = download_file(url, dest, overwrite)

    return results


def parse_ids_from_file(path: str) -> list[str]:
    ids = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Support comma- or whitespace-separated IDs on the same line
            ids.extend(token for token in line.replace(",", " ").split() if token)
    return ids


def main():
    parser = argparse.ArgumentParser(
        description="Download PDB/CIF/FASTA files from RCSB.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Format flags
    fmt = parser.add_argument_group("file formats (default: --pdb if none specified)")
    fmt.add_argument("--pdb", action="store_true", help="Download PDB format (.pdb)")
    fmt.add_argument("--cif", action="store_true", help="Download mmCIF format (.cif)")
    fmt.add_argument("--seq", action="store_true", help="Download FASTA sequence (.fasta)")

    # Input
    inp = parser.add_argument_group("input")
    inp.add_argument("ids", nargs="*", metavar="PDB_ID", help="One or more 4-character PDB IDs")
    inp.add_argument(
        "--from-file", metavar="FILE",
        help="Read PDB IDs from a file (one per line, # comments ok, commas allowed)",
    )

    # Output
    out = parser.add_argument_group("output")
    out.add_argument(
        "--outdir", metavar="DIR", default="Files",
        help="Directory to save files (default: Files)",
    )
    out.add_argument(
        "--overwrite", action="store_true",
        help="Re-download and overwrite existing files",
    )

    args = parser.parse_args()

    # Collect IDs
    ids = list(args.ids)
    if args.from_file:
        try:
            ids.extend(parse_ids_from_file(args.from_file))
        except FileNotFoundError:
            print(f"Error: file not found: {args.from_file}", file=sys.stderr)
            sys.exit(1)

    if not ids:
        parser.print_help()
        sys.exit(1)

    # Default format: PDB only
    fmt_pdb = args.pdb
    fmt_cif = args.cif
    fmt_seq = args.seq
    if not (fmt_pdb or fmt_cif or fmt_seq):
        fmt_pdb = True

    # Prepare output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Download
    ok, fail = 0, 0
    for pdb_id in ids:
        # Strip trailing _# suffix before validation (e.g. "1ABC_2" -> "1ABC")
        pdb_id = pdb_id.split("_")[0]
        if len(pdb_id) != 4 or not pdb_id.isalnum():
            print(f"\n[warn] '{pdb_id}' doesn't look like a valid PDB ID (expect 4 alphanumeric chars), skipping")
            fail += 1
            continue

        result = download_entry(pdb_id, outdir, fmt_pdb, fmt_cif, fmt_seq, args.overwrite)

        statuses = [v for v in (result["pdb"], result["cif"], result["seq"]) if v is not None]
        if all(statuses):
            ok += 1
        else:
            fail += 1

    # Summary
    total = ok + fail
    print(f"\nDone: {ok}/{total} entries downloaded successfully to '{outdir}'")
    if fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
