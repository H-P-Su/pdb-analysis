# Examples

All examples assume the virtual environment is active:

```bash
source venv/bin/activate
```

Downloaded files go to `Files/` by default. Analysis outputs go wherever you specify (examples below use the current directory or a named subdirectory).

---

## 1. Downloading Structures — `download_pdb.py`

### Download a single structure (PDB format)

```bash
python3 download_pdb.py 4HHB
```

Saves `Files/4HHB.pdb`. If the file already exists it is skipped unless `--overwrite` is passed.

### Download multiple formats at once

```bash
python3 download_pdb.py --pdb --cif --seq 4HHB 2HHB
```

Saves `4HHB.pdb`, `4HHB.cif`, `4HHB.fasta`, `2HHB.pdb`, `2HHB.cif`, `2HHB.fasta` into `Files/`.

### Download sequence only

```bash
python3 download_pdb.py --seq 4HHB 6D1Y
```

Saves `Files/4HHB.fasta` and `Files/6D1Y.fasta` (FASTA from RCSB).

### Batch download from a file

```bash
# ids.txt — one ID per line, # comments and commas are fine
python3 download_pdb.py --pdb --from-file ids.txt
```

Reads every ID in `ids.txt` and downloads PDB files for each. IDs with `_#` suffixes (e.g. `1TIM_2`) are stripped to `1TIM` automatically.

### Custom output directory

```bash
python3 download_pdb.py --pdb --outdir ./my_structures 1TIM 3J3Q
```

Saves to `my_structures/` (created if it does not exist) instead of `Files/`.

### Force re-download

```bash
python3 download_pdb.py --pdb --overwrite 4HHB
```

Re-downloads even if `4HHB.pdb` already exists.

---

## 2. Structure Summaries — `summarize_structures.py`

### Print a full metadata summary

```bash
python3 summarize_structures.py Files/4HHB.pdb
```

Prints chain count, residue count, atom count, ligands, waters, resolution, R-work/R-free, experimental method, organism, and deposition date. Example output:

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

Summarizes every PDB file in `Files/` and writes one row per structure to `summary.csv`. Columns: PDB ID, method, resolution, R-work, R-free, chains, residues, atoms, ligands, waters, organism, date.

### Missing residues report

```bash
python3 summarize_structures.py Files/6D1Y.pdb --missing
```

Compares SEQRES records (expected sequence) against ATOM records (observed residues) per chain. Reports how many residues are modeled and how many are missing.

```
Chain A: 300/320 modeled (93.8%, 20 missing)
```

### FASTA sequence from coordinates

```bash
python3 summarize_structures.py Files/4HHB.pdb --fasta
python3 summarize_structures.py Files/4HHB.pdb --fasta --chain A
```

Extracts the one-letter sequence directly from ATOM records (not SEQRES). Use `--chain` to restrict to one chain. Prints FASTA to stdout.

### Radius of gyration

```bash
python3 summarize_structures.py Files/4HHB.pdb --rg
python3 summarize_structures.py Files/4HHB.pdb --rg --chain A
```

Computes the Cα radius of gyration in Å, a measure of overall compactness. `4HHB` gives ~23.7 Å for all chains; chain A alone gives a smaller value.

### B-factor profile

```bash
python3 summarize_structures.py Files/6D1Y.pdb --bfactor
python3 summarize_structures.py Files/6D1Y.pdb --bfactor --plot --outdir ./plots
```

Prints a per-residue Cα B-factor table. Residues more than ±2σ from the chain mean are flagged as `high` (mobile/disordered) or `low` (rigid). `--plot` saves a bar chart PNG to `--outdir`.

### Ramachandran plot

```bash
python3 summarize_structures.py Files/6D1Y.pdb --ramachandran
python3 summarize_structures.py Files/6D1Y.pdb --ramachandran --plot --chain A
```

Computes φ/ψ backbone dihedral angles for every residue and classifies them as helix, sheet, allowed, or outlier. `--plot` saves a scatter plot PNG with shaded favored regions and annotated outliers. A well-determined crystal structure typically has < 5% outliers.

---

## 3. Ligand Analysis — `analyze_ligands.py`

All commands accept [PyMOL-style selections](#selection-syntax) for flexible atom group specification.

### List all ligands

```bash
python3 analyze_ligands.py Files/4HHB.pdb --ligands
python3 analyze_ligands.py Files/6D1Y.pdb --ligands
```

Lists every non-water HETATM residue: residue name, chain, residue number, and atom count. For `4HHB` this shows the four HEM groups and two phosphate ions.

### Contacts within a distance cutoff

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --contacts "resn FQJ" "polymer" --cutoff 4.5
python3 analyze_ligands.py Files/6D1Y.pdb --contacts "resn FQJ" "polymer" --out contacts.csv
```

Finds all heavy-atom pairs between the two selections within the cutoff (default 4.5 Å). Prints a table of atom1, atom2, distance. `--out` saves to CSV. For FQJ in 6D1Y this finds ~135 contacts.

### Hydrogen bonds

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --hbonds "resn FQJ" "polymer"
python3 analyze_ligands.py Files/6D1Y.pdb --hbonds "polymer" "resn FQJ" --out hbonds.csv
```

Identifies H-bonds between N/O/S donor and N/O/S acceptor heavy atoms (distance ≤ 3.5 Å, angle filter applied; no explicit hydrogens required). Run in both directions to catch all bonds. For FQJ in 6D1Y this finds ~8 H-bonds total.

### Pi interactions (stacking + cation-pi)

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --pi "polymer" "resn FQJ"
```

Detects face-to-face and edge-to-face aromatic ring stacking, and cation-pi interactions between ARG/LYS/HIS cations and aromatic rings. Protein rings are detected from known residue types (PHE, TYR, TRP, HIS); ligand rings are detected from CONECT records. For FQJ in 6D1Y this finds 5 pi interactions.

### Salt bridges

```bash
python3 analyze_ligands.py Files/4HHB.pdb --salt-bridges "polymer" "polymer"
python3 analyze_ligands.py Files/6D1Y.pdb --salt-bridges "polymer" "resn FQJ" --out sb.csv
```

Finds charged group pairs: ARG/LYS/HIS (cation) ↔ ASP/GLU (anion) within 4.0 Å. Works for protein–protein or protein–ligand. 4HHB has ~55 intra-protein salt bridges.

### Full ligand report (all interaction types)

```bash
python3 analyze_ligands.py Files/6D1Y.pdb --analyze "resn FQJ" --out fqj_report
python3 analyze_ligands.py Files/4HHB.pdb --analyze "resn HEM and chain A"
```

Runs contacts, H-bonds, pi interactions, and salt bridges in one command. Prints a summary line and, with `--out`, saves four CSVs: `fqj_report_contacts.csv`, `fqj_report_hbonds.csv`, `fqj_report_pi.csv`, `fqj_report_salt.csv`.

#### Selection syntax

Selections use PyMOL-compatible syntax:

| Selection | Meaning |
|---|---|
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

Generates two files per ligand: an interactive 3D HTML viewer and a 2D interaction fingerprint PNG.

### Interactive 3D viewer + fingerprint for any ligand

```bash
python3 visualize_interactions.py Files/6D1Y.pdb --analyze "resn FQJ"
python3 visualize_interactions.py Files/4HHB.pdb --analyze "resn HEM and chain A"
```

Saves `FQJ_A801.html` and `FQJ_A801_fingerprint.png` (named by ligand, chain, residue number) in the current directory. Open the HTML in any browser — no server or plugin needed.

**Color scheme:**
- Ligand — yellow carbons
- Contacts — cyan sticks
- H-bonds — orange sticks + yellow dashed lines
- Pi interactions — magenta sticks + magenta dashed lines
- Salt bridges — blue/red sticks + orange dashed lines

### Custom output directory and size

```bash
python3 visualize_interactions.py Files/6D1Y.pdb --analyze "resn FQJ" \
    --outdir ./images --width 1200 --height 900
```

Saves outputs to `images/` and sets the HTML viewer dimensions to 1200×900 px.

### All organic ligands at once

```bash
python3 visualize_interactions.py Files/4HHB.pdb --analyze "organic"
```

Runs the full analysis for every non-water HETATM residue in the file, generating one HTML + one PNG per ligand.

### Custom partner selection and contact cutoff

```bash
python3 visualize_interactions.py Files/6D1Y.pdb --analyze "resn FQJ" \
    --protein "polymer and chain A" --cutoff 5.0
```

Restricts the binding partner to chain A and widens the contact shell to 5 Å.

---

## 5. RMSD Analysis — `analyze_rmsd.py`

### RMSD between two structures

```bash
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A
```

Aligns sequences, superposes Cα atoms with the Kabsch algorithm, and prints the overall Cα RMSD in Å. `4HHB` vs `2HHB` (oxy vs deoxy hemoglobin, chain A) gives ~0.14 Å — nearly identical backbone.

### Per-residue RMSD

```bash
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A --per-residue
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A --per-residue \
    --plot --out per_residue.csv
```

After superposition, computes the Cα distance for each aligned residue pair. Prints the table and highlights residues above `--threshold` (default 2.0 Å). `--plot` saves a bar chart PNG; `--out` saves the table to CSV.

### Superpose two structures

```bash
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --superpose --out superposed.pdb
python3 analyze_rmsd.py Files/4HHB.pdb Files/2HHB.pdb --chain A --superpose --out superposed.pdb
```

Superposes the second structure onto the first using the Kabsch algorithm, then saves the rotated/translated coordinates to `superposed.pdb`. Prints the RMSD after superposition. The output PDB can be loaded into PyMOL or any viewer alongside the reference.

### Pairwise RMSD matrix

```bash
python3 analyze_rmsd.py Files/6D1Y.pdb Files/6D1Z.pdb Files/6D20.pdb --matrix
python3 analyze_rmsd.py Files/*.pdb --matrix --out rmsd_matrix.csv --plot
```

Computes all-vs-all pairwise Cα RMSD for N structures, printing a symmetric matrix. `--plot` saves a color-coded heatmap PNG; `--out` saves the matrix to CSV. Useful for clustering or comparing a series of related structures (e.g. a kinase in different ligand-bound states).

### NMR ensemble RMSD

```bash
python3 analyze_rmsd.py Files/1D3Z.pdb --ensemble
python3 analyze_rmsd.py Files/1D3Z.pdb --ensemble --plot --out ensemble.csv
```

For a multi-MODEL PDB file (NMR ensemble), computes the Cα RMSD of each model against the reference model (default: model 0). `1D3Z` (ubiquitin, 10 models) gives mean RMSD ~0.5 Å across the ensemble. `--plot` saves a bar chart; `--out` saves the per-model table to CSV.

### Alignment method

By default structures are matched by global sequence alignment (handles insertions and deletions). Use `--align resnum` when both structures share the same residue numbering scheme:

```bash
python3 analyze_rmsd.py Files/6D1Z.pdb Files/6D20.pdb --align resnum
```

---

## 6. Sequence Conservation — `conservation.py`

Maps per-residue conservation scores onto a PDB structure using a multi-FASTA file of homologous sequences. Uses BLOSUM62 pairwise alignment — no BLAST or network access required.

### Basic usage

```bash
python3 conservation.py Files/6D1Y.pdb homologs.fasta
```

Extracts the reference sequence from chain A ATOM records, aligns every sequence in `homologs.fasta`, and writes three output files to the same directory as the PDB:

```
6D1Y_conservation.png   # bar chart: red (variable) → blue (conserved)
6D1Y_conservation.pdb   # B-factor column = conservation score × 100
6D1Y_conservation.csv   # resnum, aa, score, entropy, gap_fraction, n_sequences
```

### Specify chain and output directory

```bash
python3 conservation.py Files/4HHB.pdb homologs.fasta --chain A --out results/
```

Uses only chain A as the reference sequence and writes outputs to `results/`.

### Skip the PDB or plot output

```bash
python3 conservation.py Files/6D1Y.pdb homologs.fasta --no-pdb     # CSV + plot only
python3 conservation.py Files/6D1Y.pdb homologs.fasta --no-plot    # CSV + PDB only
```

### Colour by conservation in PyMOL

```bash
# The script prints this command after each run:
load 6D1Y_conservation.pdb;  spectrum b, red_white_blue, minimum=0, maximum=100
```

Residues coloured white are fully variable; blue = fully conserved.

### FASTA file format

Standard multi-FASTA — one or more sequences, any number of sequences:

```
>UniProt_P12345 Human kinase
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQNTAH...
>UniProt_Q98765 Mouse kinase
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQNTSH...
```

---

## 7. Running All Examples at Once

The `examples.py` script downloads six demo structures and exercises every feature above, saving all outputs to `demo_output/`:

```bash
python3 examples.py
```

**Demo structures used:**

| ID | Description |
|---|---|
| 4HHB | Oxy-hemoglobin — HEM ligands, salt bridges, Ramachandran |
| 2HHB | Deoxy-hemoglobin — RMSD / superposition pair for 4HHB |
| 6D1Y | Kinase with inhibitor FQJ — contacts, H-bonds, pi, visualization |
| 6D1Z | Same kinase, different conformation — pairwise RMSD |
| 6D20 | Same kinase series — pairwise RMSD |
| 1D3Z | Ubiquitin NMR ensemble (10 models) — ensemble RMSD |

**Output files produced:**

| File | What it is |
|---|---|
| `batch_summary.csv` | Metadata table for all 6 structures |
| `6D1Y_bfactor.png` | B-factor bar chart for 6D1Y (kinase) |
| `6D1Y_rama.png` | Ramachandran scatter plot for 6D1Y |
| `4HHB_rama.png` | Ramachandran scatter plot for 4HHB |
| `6D1Y_FQJ_contacts.csv` | All heavy-atom contacts between FQJ and polymer |
| `6D1Y_FQJ_hbonds.csv` | H-bonds between FQJ and polymer |
| `6D1Y_FQJ_pi.csv` | Pi interactions between FQJ and polymer |
| `4HHB_salt_bridges.csv` | Protein–protein salt bridges in 4HHB |
| `demo_FQJ_A801.html` | Interactive 3D viewer for FQJ in 6D1Y |
| `demo_FQJ_A801_fingerprint.png` | 2D interaction fingerprint for FQJ |
| `demo_HEM_A142.html` | Interactive 3D viewer for HEM in 4HHB |
| `demo_HEM_A142_fingerprint.png` | 2D interaction fingerprint for HEM |
| `hhb_per_residue_rmsd.csv` | Per-residue Cα RMSD, 4HHB vs 2HHB chain A |
| `hhb_per_residue_rmsd.png` | Bar chart of per-residue RMSD |
| `2HHB_superposed.pdb` | 2HHB chain A superposed onto 4HHB chain A |
| `kinase_rmsd_matrix.png` | Pairwise RMSD heatmap for 6D1Y/6D1Z/6D20 |
| `1D3Z_ensemble_rmsd.csv` | Per-model RMSD for 1D3Z NMR ensemble |
| `1D3Z_ensemble_rmsd.png` | Bar chart of ensemble RMSD |

---

## 8. Ligand RMSD — `ligand_rmsd.py`

Computes heavy-atom RMSD for a named ligand across multiple PDB structures, after aligning each structure to a reference by backbone Cα superposition.

### Centroid reference (default)

```bash
python3 ligand_rmsd.py --ligand HEM Files/2HHB.pdb Files/4HHB.pdb
```

With no `--reference`, all structures are aligned to the first file as an anchor, and the centroid (mean coordinates of the common heavy atoms) is used as the reference. Both structures are reported symmetrically around that centroid:

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

Superimposes each structure onto the reference and measures ligand deviation from it:

```
| File     | RMSD (Å) | Z-score | Atoms matched | Included |
|----------|----------|---------|---------------|----------|
| 2HHB.pdb |    0.251 |      +1 |            43 | Yes      |
| 4HHB.pdb |    0.000 |      -1 |            43 | Yes      |

Average RMSD: 0.126 Å  (over 2 structures)
```

The 0.251 Å shift reflects the small but real rearrangement of the heme group between oxy (4HHB) and deoxy (2HHB) haemoglobin.

### Glob pattern

```bash
python3 ligand_rmsd.py --ligand ATP "structures/*.pdb"
```

Expands the glob and processes all matching files. Quotes are required to prevent the shell from expanding the pattern before the script sees it.

### Restrict alignment to one chain

```bash
python3 ligand_rmsd.py --ligand HEM --chain A Files/2HHB.pdb Files/4HHB.pdb
```

Uses only chain A Cα atoms for backbone superposition and only looks for HEM in chain A. Useful when other chains have a different ligand or when alignment is noisy across chains.

### Custom z-score threshold

```bash
python3 ligand_rmsd.py --ligand HEM --zscore 1.5 Files/2HHB.pdb Files/4HHB.pdb
```

Structures with `|z| > 1.5` are flagged as outliers and excluded from the final average. The default threshold is 2.0.

### Report z-scores without excluding anything

```bash
python3 ligand_rmsd.py --ligand HEM --no-exclude Files/2HHB.pdb Files/4HHB.pdb
```

Computes and prints z-scores but includes every structure in the average regardless of deviation.

---

## 9. Molecular Dynamics Pipeline — `run_md.py`

> **Requires GROMACS** in PATH. Install first:
> ```bash
> brew install gromacs          # macOS
> sudo apt install gromacs      # Ubuntu/Debian
> conda install -c conda-forge gromacs
> module load gromacs           # HPC cluster
> ```
>
> Required Python extras (install into venv before running):
> ```bash
> pip install pdbfixer openmm   # fills missing side-chain atoms — needed for most crystal structures
> pip install MDAnalysis        # DSSP secondary structure timeline
> ```
>
> Without `pdbfixer`, crystal structures with disordered residues (missing side-chain atoms) will
> cause `pdb2gmx` to fail. Install it first.

### Full pipeline — default 10 ns

```bash
python3 run_md.py Files/6D1Y.pdb
```

Runs all seven stages in sequence and writes everything to `Files/6D1Y_md/`:

1. **prepare** — strips HETATM/ANISOU; fills missing heavy atoms with PDBFixer
2. **setup** — builds AMBER99SB-ILDN topology, wraps in a dodecahedral TIP3P box, adds 0.15 M NaCl
3. **minimize** — steepest-descent energy minimisation to relieve steric clashes
4. **equil** — 100 ps NVT (heat to 300 K) then 100 ps NPT (equilibrate density)
5. **run** — 10 ns production MD, saving XTC every 10 ps
6. **analyze** — RMSD, RMSF, Rg, energy, H-bond count, DSSP
7. **visualize** — matplotlib plots + PyMOL session script

### Longer production run

```bash
python3 run_md.py Files/4HHB.pdb --ns 100
```

Runs a 100 ns production simulation. All other defaults remain: 300 K, 0.15 M NaCl, AMBER99SB-ILDN/TIP3P, 10 ps save interval.

### Different force field and water model

```bash
python3 run_md.py Files/4HHB.pdb --ff charmm36m-iua --water tip4p
```

CHARMM36m with TIP4P water — recommended for intrinsically disordered proteins or when accurate water diffusion matters.

### Parallel CPU run

```bash
python3 run_md.py Files/6D1Y.pdb --ncores 8
```

Pins GROMACS to 8 OpenMP threads on one MPI rank (`-ntmpi 1 -ntomp 8`). Omit `--ncores` to let GROMACS auto-detect the best thread count.

### GPU-accelerated run

```bash
python3 run_md.py Files/6D1Y.pdb --ncores 4 --gpu
```

Offloads non-bonded force calculations to GPU 0. Requires a GROMACS build with CUDA or OpenCL support.

### Run only specific stages

```bash
# Build inputs only (no simulation yet)
python3 run_md.py Files/6D1Y.pdb --steps prepare setup

# Run simulation assuming inputs are already built
python3 run_md.py Files/6D1Y.pdb --steps minimize equil run

# Post-process an existing trajectory (no re-running MD)
python3 run_md.py Files/6D1Y.pdb --steps analyze visualize
```

Every stage checks whether its output already exists and skips it if so — re-running the full command safely resumes from the last completed step.

### Custom temperature and salt

```bash
python3 run_md.py Files/4HHB.pdb --temp 310 --conc 0.05 --ns 50
```

Simulates at 310 K (physiological) with 50 mM NaCl for 50 ns.

### Fine-grained save interval

```bash
python3 run_md.py Files/6D1Y.pdb --save-ps 2
```

Saves a frame every 2 ps instead of every 10 ps — useful when studying fast conformational changes but produces larger trajectory files.

### Output files (written to `<stem>_md/`)

| File | What it is |
|---|---|
| `topol.top` | GROMACS topology (force-field parameters for the whole system) |
| `em.gro` | Energy-minimised coordinates |
| `md.tpr` | Production binary (topology + coordinates; needed for analysis) |
| `md.xtc` | Raw compressed trajectory (all atoms) |
| `md_center.xtc` | Trajectory with protein centred and PBC repaired |
| `md_fit.xtc` | Trajectory fitted to backbone (rotation/translation removed) |
| `rmsd.xvg` | Backbone RMSD vs time (nm; parse with `parse_xvg()`) |
| `rmsf.xvg` | Per-residue Cα RMSF (nm vs residue number) |
| `gyrate.xvg` | Radius of gyration vs time (nm) |
| `energy.xvg` | Potential energy, temperature, pressure vs time |
| `hbond.xvg` | Intra-protein H-bond count vs time |
| `dssp.npz` | NumPy array: secondary structure codes (frames × residues) |
| `rmsf_bfactor.pdb` | Reference PDB with RMSF (Å) in the B-factor column |
| `plots/<stem>_summary.png` | All analysis panels in one figure |
| `plots/<stem>_rmsd.png` | Backbone RMSD time series |
| `plots/<stem>_rmsf.png` | Per-residue flexibility bar chart |
| `plots/<stem>_rg.png` | Radius of gyration time series |
| `plots/<stem>_energy.png` | Energy + temperature time series |
| `plots/<stem>_hbonds.png` | H-bond count time series |
| `plots/<stem>_dssp.png` | Secondary structure content stacked area chart |
| `<stem>_session.pml` | PyMOL script — load and run with `pymol <stem>_session.pml` |
| `<stem>.pse` | PyMOL binary session saved by the .pml |

### PyMOL session

The generated `.pml` script, when opened in PyMOL, will:

1. Load the reference structure (`md_ref.pdb`) and append the fitted trajectory (`md_fit.xtc`) as additional states on the same object
2. Show cartoon representation coloured by secondary structure (helix=red, strand=yellow, loop=cyan)
3. Set up the movie with `mset 1 -N` so all N trajectory frames play in sequence
4. Save a `.pse` session file for later interactive use

```bash
# Open the session interactively
cd Files/6D1Y_md && pymol 6D1Y_session.pml

# In PyMOL — type in the command line at the bottom:
play                 # play trajectory (press Escape to stop)
set movie_fps, 10    # slow down
set movie_fps, 50    # speed up
forward              # step one frame
backward             # step back one frame
```

> **Note:** `mset 1 xN` (wrong) repeats state 1 N times — no visible motion.
> `mset 1 -N` (correct) plays states 1 through N in sequence.
> The script uses the correct form.

To render an MP4 movie (requires ffmpeg):
```bash
# In PyMOL:
import os; os.makedirs("frames", exist_ok=True)
mpng frames/frame
quit

# Then in shell:
ffmpeg -r 25 -i frames/frame%04d.png -c:v libx264 -crf 20 6D1Y_traj.mp4
```
