# PDB Toolkit

A collection of Python tools for downloading, analysing, visualising, and simulating protein structures from [RCSB PDB](https://www.rcsb.org).

| Script / App | Purpose |
|---|---|
| `app.py` | **Streamlit web app** — upload a PDB and get all analyses in a browser |
| `download_pdb.py` | Download PDB / mmCIF / FASTA files from RCSB |
| `summarize_structures.py` | Metadata, B-factors, Ramachandran, FASTA, Rg, buried surface area |
| `analyze_ligands.py` | Ligand–protein contacts, H-bonds, pi interactions, salt bridges |
| `visualize_interactions.py` | Interactive 3D HTML viewers and 2D LIGPLOT-style diagrams |
| `analyze_rmsd.py` | Cα RMSD, structure superposition, pairwise RMSD matrices |
| `ligand_rmsd.py` | Track ligand conformational changes across multiple structures |
| `conservation.py` | Map sequence conservation scores onto a PDB structure |
| `run_md.py` | Full GROMACS molecular dynamics pipeline |

---

## Quick start

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Python 3.9+ required. See [INSTALL.md](INSTALL.md) for full setup instructions including optional GROMACS for `run_md.py`.

### Launch the web app

```bash
streamlit run app.py
```

Open `http://localhost:8501`, upload any `.pdb` file, and all analyses run automatically.

---

## Web App — `app.py`

The Streamlit app provides a point-and-click interface to the three main analysis libraries. Upload a `.pdb` file via the sidebar; results appear in five tabs:

| Tab | Contents |
|-----|----------|
| **📋 Structure Summary** | Metadata table, B-factor profile, Ramachandran plot, FASTA sequences, missing residues |
| **⬛ Buried Surface Area** | Pairwise chain BSA heatmap, ranked interface bar chart, per-chain SASA table |
| **🔬 Ligand Analysis** | Contacts, H-bonds, pi interactions, salt bridges per ligand — all downloadable as CSV |
| **🌐 3D Viewer** | Embedded 3Dmol.js viewer with toggleable interaction layers; download standalone HTML |
| **🖼 2D Diagrams** | Interaction fingerprint dot-plot and LIGPLOT-style 2D diagram per ligand |

Sidebar controls: contact cutoff (Å), BSA probe radius, BSA sphere-point accuracy. All analyses re-run automatically when settings change.

---

## 1. Downloading Structures — `download_pdb.py`

No third-party dependencies — Python stdlib only.

```bash
# Single structure, PDB format (default)
python3 download_pdb.py 4HHB

# Multiple structures, all formats
python3 download_pdb.py --pdb --cif --seq 4HHB 2HHB

# Batch download from a file
python3 download_pdb.py --pdb --from-file ids.txt

# Custom output directory, force re-download
python3 download_pdb.py --pdb --outdir ./my_structures --overwrite 1TIM
```

IDs with `_#` suffixes (e.g. `1TIM_2`) are stripped to `1TIM` automatically. Files already present are skipped unless `--overwrite` is given.

### Options

| Flag | Description |
|------|-------------|
| `--pdb` | Download `.pdb` format (default if no format specified) |
| `--cif` | Download `.cif` (mmCIF) format |
| `--seq` | Download `.fasta` sequence file |
| `--from-file FILE` | Read PDB IDs from a text file (one per line; `#` comments and commas OK) |
| `--outdir DIR` | Output directory (default: `Files/`) |
| `--overwrite` | Re-download even if the file already exists |

---

## 2. Structure Summaries — `summarize_structures.py`

### Metadata summary

```bash
python3 summarize_structures.py Files/4HHB.pdb
```

```
───────────────────────────────────────────────────────
  PDB:        4HHB
  Title:      THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 ANGS
  Method:     X-RAY DIFFRACTION
  Resolution: 1.74 Å
  R-work:     0.135   R-free: None
  Organism:   Homo Sapiens
  Chains:     A,B,C,D  (A:141; B:146; C:141; D:146)
  Residues:   574   Atoms: 4779
  Ligands:    6   (HEM/A,HEM/B,HEM/C,HEM/D,PO4/B,PO4/D)
  Waters:     221
───────────────────────────────────────────────────────
```

### Batch summary to CSV

```bash
python3 summarize_structures.py Files/*.pdb --out summary.csv
```

### Missing residues

```bash
python3 summarize_structures.py Files/6D1Y.pdb --missing
```

Compares SEQRES (expected) against ATOM records (observed) per chain and reports gaps.

### FASTA from coordinates

```bash
python3 summarize_structures.py Files/4HHB.pdb --fasta
python3 summarize_structures.py Files/4HHB.pdb --fasta --chain A
```

### Radius of gyration

```bash
python3 summarize_structures.py Files/4HHB.pdb --rg
# 4HHB → Rg = 23.67 Å
```

### B-factor profile

```bash
python3 summarize_structures.py Files/6D1Y.pdb --bfactor --plot --outdir ./plots
```

Per-residue Cα B-factor with ±2σ flagging and optional bar chart PNG.

### Ramachandran plot

```bash
python3 summarize_structures.py Files/6D1Y.pdb --ramachandran --plot --chain A
```

φ/ψ backbone dihedral angles classified as helix, sheet, allowed, or outlier. Well-determined crystal structures typically have < 5% outliers.

### Buried surface area (BSA)

Importable function — not a CLI flag. Use from Python or via `app.py`:

```python
from summarize_structures import load_structure, buried_surface_areas, plot_bsa_matrix

structure = load_structure("Files/4HHB.pdb")
results = buried_surface_areas(structure, probe_radius=1.4, n_points=100)
# results sorted by bsa_total descending:
# [{'chain_1':'A','chain_2':'B','bsa_total':997.0,'interface_pct_1':13.0,...}, ...]

chain_ids = [ch.id for ch in structure[0]]
plot_bsa_matrix(results, chain_ids, "4HHB_bsa.png", title="4HHB")
```

`BSA(A, B) = ( SASA_A + SASA_B − SASA_AB ) / 2` using the Shrake-Rupley rolling-sphere algorithm. Each result dict contains:

| Field | Description |
|-------|-------------|
| `chain_1`, `chain_2` | Chain IDs |
| `sasa_1`, `sasa_2` | Isolated chain SASA (Å²) |
| `sasa_complex` | Combined SASA when both chains present (Å²) |
| `bsa_total` | Total buried interface area (Å²) |
| `bsa_on_1`, `bsa_on_2` | Surface buried on each individual chain (Å²) |
| `interface_pct_1`, `interface_pct_2` | % of each chain's SASA buried at this interface |

**4HHB example results** (human deoxyhaemoglobin α₂β₂ tetramer):

| Chains | BSA (Å²) | Notes |
|--------|----------|-------|
| A–B | 997 | α₁β₁ dimer interface |
| C–D | 939 | α₂β₂ dimer interface |
| A–D | 806 | α₁β₂ interface (weaker) |
| B–C | 805 | α₂β₁ interface (weaker) |
| A–C | 322 | Minimal contact |
| B–D | 0 | No contact |

---

## 3. Ligand Analysis — `analyze_ligands.py`

### List all ligands

```bash
python3 analyze_ligands.py Files/4HHB.pdb --ligands
python3 analyze_ligands.py Files/6D1Y.pdb --ligands
```

### Contacts, H-bonds, pi interactions, salt bridges

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --contacts "resn FQJ" "polymer" --cutoff 4.5
python3 analyze_ligands.py Files/6D1Y.pdb --hbonds "resn FQJ" "polymer"
python3 analyze_ligands.py Files/6D1Y.pdb --pi "polymer" "resn FQJ"
python3 analyze_ligands.py Files/4HHB.pdb --salt-bridges "polymer" "polymer"
```

### Full ligand report — auto-detect primary ligands (default)

```bash
python3 analyze_ligands.py Files/6D1Y.pdb
python3 analyze_ligands.py Files/6D1Y.pdb --out fqj_report
```

When no mode flag is given, `primary_ligands()` filters out water, excipients (SO4, GOL, EDO, ions, buffers, …), and residues with fewer than 7 heavy atoms, then runs a full analysis on what remains. With `--out`, four CSVs are saved per ligand: `_contacts.csv`, `_hbonds.csv`, `_pi.csv`, `_salt.csv`.

### Selection syntax

| Selection | Meaning |
|-----------|---------|
| `polymer` | All protein and nucleic acid residues |
| `hetatm` | All HETATM records |
| `organic` | HETATM excluding water |
| `resn ATP` | Residue named ATP |
| `resn PHE+TYR+TRP` | Any of these names |
| `chain A+B` | Chains A and B |
| `resi 100-200` | Residue numbers 100–200 |
| `b > 50` | Atoms with B-factor > 50 |
| `within 5 of resn ATP` | Anything within 5 Å of ATP |
| `not solvent` | Exclude water |

---

## 4. Visualization — `visualize_interactions.py`

Generates five output files per ligand:

| File | Contents |
|------|----------|
| `{LIG}_{chain}{resi}.html` | Interactive 3D viewer — all interaction types toggleable |
| `{LIG}_{chain}{resi}_fingerprint.png` | 2D dot-plot (residues × interaction types) |
| `{LIG}_{chain}{resi}_2d.png` | LIGPLOT-style 2D diagram |
| `{LIG}_{chain}{resi}_{type}.html` | Per-type 3D viewer (one per detected type) |
| `{LIG}_{chain}{resi}_{type}_fingerprint.png` | Per-type filtered fingerprint |

```bash
# Auto-detect primary ligands
python3 visualize_interactions.py Files/6D1Y.pdb

# Explicit ligand selection
python3 visualize_interactions.py Files/6D1Y.pdb --analyze "resn FQJ"
python3 visualize_interactions.py Files/4HHB.pdb --analyze "resn HEM and chain A"

# Custom output directory and viewer size
python3 visualize_interactions.py Files/6D1Y.pdb --analyze "resn FQJ" \
    --outdir ./images --width 1200 --height 900
```

Open any `.html` file in a browser — no server, no plugin, no py3Dmol needed. The file is self-contained and loads 3Dmol.js from CDN.

### Color scheme

| Element | Color |
|---------|-------|
| Ligand | Yellow carbons |
| Contacts | Cyan sticks; spoked arcs in 2D |
| H-bonds | Orange sticks + yellow dashed lines |
| Pi interactions | Magenta sticks + magenta dashed lines |
| Salt bridges | Blue (+) / red (−) sticks + orange dashed lines |

---

## 5. RMSD Analysis — `analyze_rmsd.py`

```bash
# Overall Cα RMSD
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A
# → 4HHB vs 2HHB chain A ≈ 0.14 Å

# Per-residue RMSD with plot
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A --per-residue --plot

# Superpose and save
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --superpose --out superposed.pdb

# All-vs-all pairwise matrix
python3 analyze_rmsd.py Files/*.pdb --matrix --out rmsd_matrix.csv --plot

# NMR ensemble
python3 analyze_rmsd.py Files/1D3Z.pdb --ensemble --plot
# → 1D3Z (ubiquitin, 10 models): mean RMSD ≈ 0.5 Å
```

By default, structures are matched by global sequence alignment. Use `--align resnum` when residue numbering is shared.

---

## 6. Ligand RMSD — `ligand_rmsd.py`

Tracks heavy-atom RMSD for a named ligand across multiple structures after backbone Cα superposition.

```bash
# Centroid reference (default)
python3 ligand_rmsd.py --ligand HEM Files/2HHB.pdb Files/4HHB.pdb

# Explicit reference
python3 ligand_rmsd.py --ligand HEM --reference Files/4HHB.pdb Files/2HHB.pdb Files/4HHB.pdb
# → 2HHB deviates 0.251 Å from 4HHB — real oxy/deoxy heme rearrangement

# Glob pattern
python3 ligand_rmsd.py --ligand ATP "structures/*.pdb"

# Restrict to one chain, custom z-score threshold
python3 ligand_rmsd.py --ligand HEM --chain A --zscore 1.5 Files/2HHB.pdb Files/4HHB.pdb
```

---

## 7. Sequence Conservation — `conservation.py`

Maps per-residue conservation scores onto a PDB structure from a multi-FASTA alignment. Uses BLOSUM62 — no BLAST or network access required.

```bash
python3 conservation.py Files/6D1Y.pdb homologs.fasta
python3 conservation.py Files/4HHB.pdb homologs.fasta --chain A --out results/
```

**Outputs:** `<stem>_conservation.png` (bar chart, red→blue), `<stem>_conservation.pdb` (scores in B-factor column), `<stem>_conservation.csv`.

**Colour by conservation in PyMOL:**
```
load 6D1Y_conservation.pdb;  spectrum b, red_white_blue, minimum=0, maximum=100
```

---

## 8. Molecular Dynamics Pipeline — `run_md.py`

> **Requires GROMACS** in PATH — see [INSTALL.md](INSTALL.md).

```bash
# Full 10 ns pipeline
python3 run_md.py Files/6D1Y.pdb

# 100 ns run
python3 run_md.py Files/4HHB.pdb --ns 100

# CHARMM36m force field
python3 run_md.py Files/4HHB.pdb --ff charmm36m-iua --water tip4p

# GPU-accelerated, 8 cores
python3 run_md.py Files/6D1Y.pdb --ncores 8 --gpu

# Build inputs only (no simulation)
python3 run_md.py Files/6D1Y.pdb --steps prepare setup

# Post-process an existing trajectory
python3 run_md.py Files/6D1Y.pdb --steps analyze visualize
```

Seven stages: **prepare → setup → minimize → equil → run → analyze → visualize**. Each stage checks for existing outputs and skips if complete — re-running safely resumes from the last completed step.

### Key outputs (`<stem>_md/`)

| File | Contents |
|------|----------|
| `md.xtc` | Raw trajectory |
| `md_fit.xtc` | Backbone-fitted trajectory |
| `rmsd.xvg` | Backbone RMSD vs time |
| `rmsf.xvg` | Per-residue Cα RMSF |
| `gyrate.xvg` | Radius of gyration vs time |
| `dssp.npz` | Secondary structure array (frames × residues) |
| `rmsf_bfactor.pdb` | PDB with RMSF in B-factor column |
| `plots/<stem>_summary.png` | All analysis panels |
| `<stem>_session.pml` | PyMOL session script |
| `<stem>_report.html` | Self-contained HTML report |

### PyMOL playback

```bash
cd Files/6D1Y_md && pymol 6D1Y_session.pml
# In PyMOL: play  (Escape to stop)
```

> **Note:** use `mset 1 -N` (not `mset 1 xN`) to play all N frames sequentially.

---

## Running all examples

```bash
python3 examples.py
```

Downloads six demo structures (4HHB, 2HHB, 6D1Y, 6D1Z, 6D20, 1D3Z) and exercises every feature, saving all outputs to `demo_output/`.

---

## Importable API

All scripts expose their core functions for use in other projects:

```python
from summarize_structures import (
    load_structure, summarize_structure, bfactor_stats, ramachandran_data,
    radius_of_gyration, extract_fasta, find_missing_residues,
    buried_surface_areas, plot_bsa_matrix,
    plot_bfactor, plot_ramachandran,
)

from analyze_ligands import (
    load_structure, get_ligands, primary_ligands,
    find_contacts, find_hydrogen_bonds, find_pi_interactions, find_salt_bridges,
    analyze_ligand,
)

from visualize_interactions import (
    save_html, plot_interaction_summary, plot_ligand_2d,
    save_type_images, visualize,
)
```
