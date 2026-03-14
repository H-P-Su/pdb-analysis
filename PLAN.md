# PDB Downloader — Project Plan

## Goal
A standalone Python CLI script to download structure and sequence files from RCSB (rcsb.org) given one or more PDB IDs.

## Project Structure
```
PDBs/
├── PLAN.md             # This file
├── download_pdb.py     # Main CLI script
├── venv/               # Python virtual environment
└── Files/              # Default download output directory
```

## Setup
```bash
python3 -m venv venv
source venv/bin/activate
```
No third-party dependencies — stdlib only (`argparse`, `urllib`, `pathlib`).

## Script: download_pdb.py

### Requirements
- Download files from RCSB using public URLs (no API key needed)
- Supported formats:
  - `--pdb` → `https://files.rcsb.org/download/{ID}.pdb`
  - `--cif` → `https://files.rcsb.org/download/{ID}.cif`
  - `--seq` → `https://www.rcsb.org/fasta/entry/{ID}` (saves as `.fasta`)
- Default format is `--pdb` if no format flag is specified
- Accept multiple PDB IDs as positional arguments
- Accept a file of IDs via `--from-file` (one per line, `#` comments, comma-separated ok)
- Default output directory: `Files/` (created automatically if missing)
- Skip existing files unless `--overwrite` is passed
- Print per-file status and a summary line at the end
- Exit with code 1 if any entry fails

### ID Handling
- Uppercase all IDs
- Strip trailing `_#` suffixes before processing (e.g. `1ABC_2` → `1ABC`)
- Validate that the result is exactly 4 alphanumeric characters; warn and skip otherwise

### CLI Interface
```
usage: download_pdb.py [--pdb] [--cif] [--seq]
                       [--from-file FILE]
                       [--outdir DIR] [--overwrite]
                       PDB_ID [PDB_ID ...]
```

### Key Functions
| Function | Purpose |
|---|---|
| `download_file(url, dest, overwrite)` | Fetch a single URL to disk; handles 404 and network errors |
| `download_entry(pdb_id, outdir, fmt_pdb, fmt_cif, fmt_seq, overwrite)` | Download all requested formats for one PDB ID |
| `parse_ids_from_file(path)` | Read IDs from a flat text file |
| `main()` | Argument parsing, ID normalization, orchestration |

## Usage Examples
```bash
# Single entry, PDB only (default)
python3 download_pdb.py 1TIM

# Multiple entries, all formats
python3 download_pdb.py --pdb --cif --seq 1TIM 4HHB

# IDs with suffixes are handled automatically
python3 download_pdb.py 1TIM_2 4HHB_10

# Batch from file
python3 download_pdb.py --cif --from-file ids.txt

# Custom output dir, force re-download
python3 download_pdb.py --pdb --outdir ./my_structs --overwrite 1TIM
```
