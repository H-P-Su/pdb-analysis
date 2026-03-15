# PDB Toolkit

A collection of eight command-line Python scripts for downloading, analysing, visualising, and simulating protein structures from [RCSB PDB](https://www.rcsb.org).

| Script | Purpose |
|--------|---------|
| `download_pdb.py` | Download PDB / mmCIF / FASTA files from RCSB |
| `summarize_structures.py` | Extract metadata, geometry, B-factors, Ramachandran analysis |
| `analyze_ligands.py` | Detect ligand–protein contacts, H-bonds, pi interactions, salt bridges |
| `visualize_interactions.py` | Generate interactive 3D HTML viewers and 2D interaction diagrams |
| `analyze_rmsd.py` | Compute Cα RMSD, superpose structures, build pairwise RMSD matrices |
| `ligand_rmsd.py` | Track ligand conformational changes across multiple structures |
| `conservation.py` | Map sequence conservation scores onto a PDB structure |
| `run_md.py` | Run a full GROMACS molecular dynamics pipeline |

---

## Quick start

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Python 3.9+ required. GROMACS must be installed separately to use `run_md.py` (see that section).

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
| `--from-file FILE` | Read PDB IDs from a text file |
| `--outdir DIR` | Output directory (default: `Files/`) |
| `--overwrite` | Re-download even if the file already exists |

### Batch file format

```
# This is a comment
1TIM
4HHB, 2XYZ
3DEF
```

---

## 2. Structure Summaries — `summarize_structures.py`

### Print a full metadata summary

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

Writes one row per structure. Columns: PDB ID, method, resolution, R-work, R-free, chains, residues, atoms, ligands, waters, organism, date.

### Missing residues report

```bash
python3 summarize_structures.py Files/6D1Y.pdb --missing
```

Compares SEQRES (expected) against ATOM records (observed) per chain:

```
Chain A: 300/320 modeled (93.8%, 20 missing)
```

### FASTA sequence from coordinates

```bash
python3 summarize_structures.py Files/4HHB.pdb --fasta
python3 summarize_structures.py Files/4HHB.pdb --fasta --chain A
```

Extracts the one-letter sequence from ATOM records (not SEQRES). Prints FASTA to stdout.

### Radius of gyration

```bash
python3 summarize_structures.py Files/4HHB.pdb --rg
python3 summarize_structures.py Files/4HHB.pdb --rg --chain A
```

Computes the Cα radius of gyration in Å. `4HHB` gives ~23.7 Å for all chains.

### B-factor profile

```bash
python3 summarize_structures.py Files/6D1Y.pdb --bfactor
python3 summarize_structures.py Files/6D1Y.pdb --bfactor --plot --outdir ./plots
```

Per-residue Cα B-factor table with residues flagged `high`/`low` when more than ±2σ from the chain mean. `--plot` saves a bar chart PNG.

### Ramachandran plot

```bash
python3 summarize_structures.py Files/6D1Y.pdb --ramachandran
python3 summarize_structures.py Files/6D1Y.pdb --ramachandran --plot --chain A
```

Computes φ/ψ backbone dihedral angles and classifies each residue as helix, sheet, allowed, or outlier. `--plot` saves a scatter plot with shaded favoured regions and annotated outliers. Well-determined crystal structures typically have < 5% outliers.

---

## 3. Ligand Analysis — `analyze_ligands.py`

### List all ligands

```bash
python3 analyze_ligands.py Files/4HHB.pdb --ligands
python3 analyze_ligands.py Files/6D1Y.pdb --ligands
```

Lists every non-water HETATM residue: residue name, chain, residue number, atom count.

### Contacts within a distance cutoff

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --contacts "resn FQJ" "polymer" --cutoff 4.5
python3 analyze_ligands.py Files/6D1Y.pdb --contacts "resn FQJ" "polymer" --out contacts.csv
```

All heavy-atom pairs between the two selections within the cutoff (default 4.5 Å). For FQJ in 6D1Y this finds ~135 contacts.

### Hydrogen bonds

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --hbonds "resn FQJ" "polymer"
python3 analyze_ligands.py Files/6D1Y.pdb --hbonds "polymer" "resn FQJ" --out hbonds.csv
```

Identifies H-bonds between N/O/S donor and N/O/S acceptor heavy atoms (distance ≤ 3.5 Å, angle filter; no explicit hydrogens required). Run in both directions to catch all bonds.

### Pi interactions (stacking + cation-pi)

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --pi "polymer" "resn FQJ"
```

Detects face-to-face and edge-to-face aromatic ring stacking, and cation-pi interactions between ARG/LYS/HIS and aromatic rings. For FQJ in 6D1Y this finds 5 pi interactions.

### Salt bridges

```bash
python3 analyze_ligands.py Files/4HHB.pdb --salt-bridges "polymer" "polymer"
python3 analyze_ligands.py Files/6D1Y.pdb --salt-bridges "polymer" "resn FQJ" --out sb.csv
```

ARG/LYS/HIS (cation) ↔ ASP/GLU (anion) pairs within 4.0 Å. Works for protein–protein or protein–ligand. 4HHB has ~55 intra-protein salt bridges.

### Full ligand report — auto-detect (default)

```bash
python3 analyze_ligands.py Files/6D1Y.pdb
python3 analyze_ligands.py Files/6D1Y.pdb --out fqj_report
```

When no mode flag is given, the script auto-detects primary ligands by filtering out water, common crystallographic excipients (SO4, GOL, EDO, buffer molecules, ions), and small residues (< 7 heavy atoms). A full analysis is run on what remains, largest first. With `--out`, four CSVs are saved: `fqj_report_contacts.csv`, `fqj_report_hbonds.csv`, `fqj_report_pi.csv`, `fqj_report_salt.csv`.

### Full ligand report — explicit selection

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --analyze "resn FQJ" --out fqj_report
python3 analyze_ligands.py Files/4HHB.pdb --analyze "resn HEM and chain A"
```

### Selection syntax

Selections use PyMOL-compatible syntax:

| Selection | Meaning |
|-----------|---------|
| `polymer` | All protein and nucleic acid residues |
| `hetatm` | All HETATM records |
| `organic` | HETATM excluding water |
| `solvent` | Water molecules only |
| `resn ATP` | Residue named ATP |
| `resn PHE+TYR+TRP` | Any of these residue names |
| `chain A` | Chain A only |
| `chain A+B` | Chains A and B |
| `resi 100` | Residue number 100 |
| `resi 100-200` | Residue numbers 100–200 |
| `name CA` | Atoms named CA |
| `elem N+O` | Nitrogen or oxygen atoms |
| `b > 50` | Atoms with B-factor above 50 |
| `within 5 of resn ATP` | Anything within 5 Å of ATP |
| `not solvent` | Exclude water |
| `resn ATP and chain B` | ATP in chain B only |

---

## 4. Visualization — `visualize_interactions.py`

Generates five output files per ligand per run:

| File | Contents |
|------|----------|
| `{LIG}_{chain}{resi}.html` | Interactive 3D viewer — all interaction types toggleable |
| `{LIG}_{chain}{resi}_fingerprint.png` | 2D dot-plot (residues × interaction types) |
| `{LIG}_{chain}{resi}_2d.png` | LIGPLOT-style 2D diagram |
| `{LIG}_{chain}{resi}_{type}.html` | Per-type 3D viewer (one per type with data) |
| `{LIG}_{chain}{resi}_{type}_fingerprint.png` | Per-type filtered fingerprint |

### Auto-detect primary ligands (default)

```bash
python3 visualize_interactions.py Files/6D1Y.pdb
python3 visualize_interactions.py Files/6D1Y.pdb --outdir ./images
```

Filters out excipients/ions/small residues and processes the remaining ligands. Prints `Auto-detected primary ligand(s): …` to confirm.

### Explicit ligand selection

```bash
python3 visualize_interactions.py Files/6D1Y.pdb --analyze "resn FQJ"
python3 visualize_interactions.py Files/4HHB.pdb --analyze "resn HEM and chain A"
python3 visualize_interactions.py Files/4HHB.pdb --analyze "organic"
```

Open any `.html` file in a browser — no server or plugin needed.

### Color scheme

| Element | Color |
|---------|-------|
| Ligand | Yellow carbons |
| Contacts | Cyan sticks; spoked arcs in 2D diagram |
| H-bonds | Orange sticks + yellow dashed lines |
| Pi interactions | Magenta sticks + magenta dashed lines |
| Salt bridges | Blue/red sticks + orange dashed lines |

### LIGPLOT-style 2D diagram (`*_2d.png`)

- Ligand atoms projected to their best-fit plane via SVD, drawn as a stick diagram with CPK colouring (C=dark gray, N=blue, O=red, S=yellow, P=orange)
- Covalent bonds inferred from 3D distance (< 1.9 Å)
- Protein residues placed evenly around a ring as labelled colour-coded boxes
- Spoked arcs for pure hydrophobic contacts; dashed lines for H-bonds, pi, and salt bridges

### Custom viewer size and partner selection

```bash
python3 visualize_interactions.py Files/6D1Y.pdb --analyze "resn FQJ" \
    --outdir ./images --width 1200 --height 900

python3 visualize_interactions.py Files/6D1Y.pdb --analyze "resn FQJ" \
    --protein "polymer and chain A" --cutoff 5.0
```

---

## 5. RMSD Analysis — `analyze_rmsd.py`

### RMSD between two structures

```bash
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A
```

Aligns sequences, superposes Cα atoms (Kabsch algorithm), and prints overall Cα RMSD. `4HHB` vs `2HHB` chain A gives ~0.14 Å.

### Per-residue RMSD

```bash
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A --per-residue
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A --per-residue \
    --plot --out per_residue.csv
```

Cα distance for each aligned residue pair after superposition. Residues above `--threshold` (default 2.0 Å) are highlighted. `--plot` saves a bar chart PNG.

### Superpose two structures

```bash
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --superpose --out superposed.pdb
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A --superpose --out superposed.pdb
```

Rotates/translates the second structure onto the first and saves the result to `superposed.pdb`. Loadable directly in PyMOL alongside the reference.

### Pairwise RMSD matrix

```bash
python3 analyze_rmsd.py Files/6D1Y.pdb Files/6D1Z.pdb Files/6D20.pdb --matrix
python3 analyze_rmsd.py Files/*.pdb --matrix --out rmsd_matrix.csv --plot
```

All-vs-all pairwise Cα RMSD for N structures. `--plot` saves a colour-coded heatmap PNG.

### NMR ensemble RMSD

```bash
python3 analyze_rmsd.py Files/1D3Z.pdb --ensemble
python3 analyze_rmsd.py Files/1D3Z.pdb --ensemble --plot --out ensemble.csv
```

For a multi-MODEL PDB (NMR ensemble), computes Cα RMSD of each model against model 0. `1D3Z` (ubiquitin, 10 models) gives mean RMSD ~0.5 Å.

### Alignment method

By default, structures are matched by global sequence alignment (handles insertions/deletions). Use `--align resnum` when both structures share the same residue numbering:

```bash
python3 analyze_rmsd.py Files/6D1Z.pdb Files/6D20.pdb --align resnum
```

---

## 6. Ligand RMSD — `ligand_rmsd.py`

Tracks heavy-atom RMSD for a named ligand across multiple structures, after aligning each to a reference by backbone Cα superposition.

### Centroid reference (default)

```bash
python3 ligand_rmsd.py --ligand HEM Files/2HHB.pdb Files/4HHB.pdb
```

```
| File     | RMSD (Å) | Z-score | Atoms matched | Included |
|----------|----------|---------|---------------|----------|
| 2HHB.pdb |    0.126 |       0 |            43 | Yes      |
| 4HHB.pdb |    0.126 |       0 |            43 | Yes      |

Average RMSD: 0.126 Å  (over 2 structures)
```

### Explicit reference

```bash
python3 ligand_rmsd.py --ligand HEM --reference Files/4HHB.pdb Files/2HHB.pdb Files/4HHB.pdb
```

Measures each structure's ligand deviation from the specified reference. The 0.251 Å shift above reflects the real heme rearrangement between oxy and deoxy haemoglobin.

### Glob pattern

```bash
python3 ligand_rmsd.py --ligand ATP "structures/*.pdb"
```

Quotes are required to prevent the shell from expanding the pattern before the script sees it.

### Restrict alignment to one chain

```bash
python3 ligand_rmsd.py --ligand HEM --chain A Files/2HHB.pdb Files/4HHB.pdb
```

Uses only chain A Cα atoms for superposition and only searches for the ligand in chain A.

### Z-score filtering

```bash
# Custom threshold (default is 2.0)
python3 ligand_rmsd.py --ligand HEM --zscore 1.5 Files/2HHB.pdb Files/4HHB.pdb

# Report z-scores without excluding outliers
python3 ligand_rmsd.py --ligand HEM --no-exclude Files/2HHB.pdb Files/4HHB.pdb
```

---

## 7. Sequence Conservation — `conservation.py`

Maps per-residue conservation scores onto a PDB structure using a multi-FASTA of homologous sequences. Uses BLOSUM62 pairwise alignment — no BLAST or network access required.

### Basic usage

```bash
python3 conservation.py Files/6D1Y.pdb homologs.fasta
```

Extracts the reference sequence from ATOM records, aligns every sequence in `homologs.fasta`, and writes:

```
6D1Y_conservation.png   # bar chart: red (variable) → blue (conserved)
6D1Y_conservation.pdb   # B-factor column = conservation score × 100
6D1Y_conservation.csv   # resnum, aa, score, entropy, gap_fraction, n_sequences
```

### Specify chain and output directory

```bash
python3 conservation.py Files/4HHB.pdb homologs.fasta --chain A --out results/
```

### Skip specific outputs

```bash
python3 conservation.py Files/6D1Y.pdb homologs.fasta --no-pdb     # CSV + plot only
python3 conservation.py Files/6D1Y.pdb homologs.fasta --no-plot    # CSV + PDB only
```

### Colour by conservation in PyMOL

After running, the script prints:

```bash
load 6D1Y_conservation.pdb;  spectrum b, red_white_blue, minimum=0, maximum=100
```

White = fully variable; blue = fully conserved.

### FASTA file format

Standard multi-FASTA — any number of sequences:

```
>UniProt_P12345 Human kinase
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQNTAH...
>UniProt_Q98765 Mouse kinase
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQNTSH...
```

---

## 8. Molecular Dynamics Pipeline — `run_md.py`

> **Requires GROMACS** in PATH. Install first:
> ```bash
> brew install gromacs          # macOS
> sudo apt install gromacs      # Ubuntu/Debian
> conda install -c conda-forge gromacs
> module load gromacs           # HPC cluster
> ```
>
> Required Python extras:
> ```bash
> pip install pdbfixer openmm   # fills missing side-chain atoms
> pip install MDAnalysis        # DSSP secondary structure timeline
> ```

### Full pipeline — default 10 ns

```bash
python3 run_md.py Files/6D1Y.pdb
```

Runs all seven stages and writes everything to `Files/6D1Y_md/`:

1. **prepare** — strips HETATM/ANISOU; fills missing heavy atoms with PDBFixer
2. **setup** — builds AMBER99SB-ILDN topology, dodecahedral TIP3P box, 0.15 M NaCl
3. **minimize** — steepest-descent energy minimisation
4. **equil** — 100 ps NVT (heat to 300 K) then 100 ps NPT (equilibrate density)
5. **run** — 10 ns production MD, saving XTC every 10 ps
6. **analyze** — RMSD, RMSF, Rg, energy, H-bond count, DSSP
7. **visualize** — Matplotlib plots + PyMOL session script

Each stage checks whether its output already exists and skips it if so — re-running the full command safely resumes from the last completed step.

### Common options

```bash
# Longer run
python3 run_md.py Files/4HHB.pdb --ns 100

# Different force field and water model
python3 run_md.py Files/4HHB.pdb --ff charmm36m-iua --water tip4p

# Parallel CPU run
python3 run_md.py Files/6D1Y.pdb --ncores 8

# GPU-accelerated run
python3 run_md.py Files/6D1Y.pdb --ncores 4 --gpu

# Physiological conditions
python3 run_md.py Files/4HHB.pdb --temp 310 --conc 0.05 --ns 50

# Save frames every 2 ps (default: 10 ps)
python3 run_md.py Files/6D1Y.pdb --save-ps 2
```

### Run only specific stages

```bash
# Build inputs only (no simulation yet)
python3 run_md.py Files/6D1Y.pdb --steps prepare setup

# Run simulation assuming inputs are already built
python3 run_md.py Files/6D1Y.pdb --steps minimize equil run

# Post-process an existing trajectory (no re-running MD)
python3 run_md.py Files/6D1Y.pdb --steps analyze visualize
```

### Key output files (`<stem>_md/`)

| File | Contents |
|------|----------|
| `topol.top` | GROMACS topology |
| `em.gro` | Energy-minimised coordinates |
| `md.xtc` | Raw compressed trajectory |
| `md_fit.xtc` | Trajectory fitted to backbone |
| `rmsd.xvg` | Backbone RMSD vs time |
| `rmsf.xvg` | Per-residue Cα RMSF |
| `gyrate.xvg` | Radius of gyration vs time |
| `energy.xvg` | Potential energy, temperature, pressure |
| `hbond.xvg` | Intra-protein H-bond count vs time |
| `dssp.npz` | Secondary structure array (frames × residues) |
| `rmsf_bfactor.pdb` | Reference PDB with RMSF in B-factor column |
| `plots/<stem>_summary.png` | All analysis panels in one figure |
| `<stem>_session.pml` | PyMOL script — load with `pymol <stem>_session.pml` |
| `<stem>.pse` | PyMOL binary session |

### PyMOL session

```bash
cd Files/6D1Y_md && pymol 6D1Y_session.pml

# In PyMOL:
play                 # play trajectory (Escape to stop)
set movie_fps, 10    # slow down
set movie_fps, 50    # speed up
forward              # step one frame
backward             # step back one frame
```

To render an MP4 (requires ffmpeg):

```bash
# In PyMOL:
import os; os.makedirs("frames", exist_ok=True)
mpng frames/frame
quit

# In shell:
ffmpeg -r 25 -i frames/frame%04d.png -c:v libx264 -crf 20 6D1Y_traj.mp4
```

---

## Running all examples at once

`examples.py` downloads six demo structures and exercises every feature above, saving all outputs to `demo_output/`:

```bash
python3 examples.py
```

**Demo structures:**

| ID | Description |
|----|-------------|
| 4HHB | Oxy-haemoglobin — HEM ligands, salt bridges, Ramachandran |
| 2HHB | Deoxy-haemoglobin — RMSD / superposition pair for 4HHB |
| 6D1Y | Kinase with inhibitor FQJ — contacts, H-bonds, pi, visualization |
| 6D1Z | Same kinase, different conformation — pairwise RMSD |
| 6D20 | Same kinase series — pairwise RMSD |
| 1D3Z | Ubiquitin NMR ensemble (10 models) — ensemble RMSD |
