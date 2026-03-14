# PDB Downloader

A lightweight command-line tool to download protein structure and sequence files from [RCSB PDB](https://www.rcsb.org).

No third-party dependencies — Python stdlib only.

## Features

- Download **PDB**, **mmCIF**, and/or **FASTA sequence** files
- Process **multiple PDB IDs** in one command
- Load IDs from a **batch file** (supports comments and comma-separated entries)
- Automatically strips `_#` suffixes from IDs (e.g. `1ABC_2` → `1ABC`)
- Skips already-downloaded files unless `--overwrite` is specified
- Saves to `Files/` by default

## Usage

```bash
python3 download_pdb.py [OPTIONS] PDB_ID [PDB_ID ...]
```

### Options

| Flag | Description |
|------|-------------|
| `--pdb` | Download PDB format `.pdb` (default if no format specified) |
| `--cif` | Download mmCIF format `.cif` |
| `--seq` | Download FASTA sequence `.fasta` |
| `--from-file FILE` | Read PDB IDs from a file |
| `--outdir DIR` | Output directory (default: `Files/`) |
| `--overwrite` | Re-download and overwrite existing files |

### Examples

```bash
# Single entry, PDB format (default)
python3 download_pdb.py 1TIM

# Multiple entries, all formats
python3 download_pdb.py --pdb --cif --seq 1TIM 4HHB

# IDs with suffixes are handled automatically
python3 download_pdb.py 1TIM_2 4HHB_10

# CIF only, from a batch file
python3 download_pdb.py --cif --from-file ids.txt

# Custom output directory, force re-download
python3 download_pdb.py --pdb --outdir ./my_structs --overwrite 1TIM
```

### Batch file format

```
# This is a comment
1TIM
4HHB, 2XYZ
3DEF
```

## File Sources

| Format | URL |
|--------|-----|
| `.pdb` | `https://files.rcsb.org/download/{ID}.pdb` |
| `.cif` | `https://files.rcsb.org/download/{ID}.cif` |
| `.fasta` | `https://www.rcsb.org/fasta/entry/{ID}` |
