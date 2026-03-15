# PDB Toolkit — Project Plan

## Goal
A suite of standalone Python CLI scripts for downloading, analysing, visualising, and running molecular dynamics on protein structures from RCSB.

## Project Structure
```
PDBs/
├── PLAN.md                    # This file
├── MD_TUTORIAL.md             # Science behind the MD pipeline
├── examples.md                # Runnable CLI examples for every script
├── future_features.md         # Feature backlog
├── requirements.txt           # Python dependencies
├── download_pdb.py            # Download PDB/CIF/FASTA from RCSB
├── summarize_structures.py    # Metadata, B-factor, Ramachandran, FASTA, Rg, batch CSV
├── analyze_ligands.py         # Ligand contacts, H-bonds, pi, salt bridges
├── visualize_interactions.py  # 3D HTML viewer + 2D fingerprint PNG
├── analyze_rmsd.py            # Cα RMSD, superposition, pairwise matrix, NMR ensemble
├── ligand_rmsd.py             # Ligand heavy-atom RMSD across multiple structures
├── conservation.py            # Per-residue sequence conservation from FASTA
├── run_md.py                  # Full GROMACS MD pipeline
├── examples.py                # Runnable demo of all scripts
├── venv/                      # Python virtual environment (not in git)
└── Files/                     # Downloaded structures and outputs (not in git)
```

## Setup
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

## Script: download_pdb.py

### Requirements
- stdlib only (`argparse`, `urllib`, `pathlib`) — no third-party dependencies
- Supported formats: `--pdb`, `--cif`, `--seq` (FASTA)
- Accept multiple IDs as positional args or `--from-file`
- Default output directory: `Files/`; skip existing unless `--overwrite`
- Strip `_#` suffixes; validate 4-character alphanumeric IDs

### CLI
```bash
python3 download_pdb.py 4HHB
python3 download_pdb.py --pdb --cif --seq 4HHB 2HHB
python3 download_pdb.py --from-file ids.txt --outdir ./my_structs
python3 download_pdb.py --overwrite 4HHB
```

---

## Script: summarize_structures.py

### Dependencies
BioPython, NumPy, Matplotlib

### Key Functions
| Function | Description |
|---|---|
| `parse_header(path)` | Resolution, R-work/free, method, organism, date |
| `summarize_structure(path)` | Full summary dict |
| `batch_summary(paths, out_csv)` | Summarize N structures → CSV |
| `find_missing_residues(path)` | SEQRES vs ATOM comparison per chain |
| `extract_fasta(path, chain)` | FASTA from ATOM coordinates |
| `bfactor_stats(structure, chain)` | Per-residue Cα B-factor with ±2σ flags |
| `radius_of_gyration(structure, chain)` | Cα Rg in Å |
| `ramachandran_data(structure, chain)` | phi/psi angles + region classification |
| `plot_bfactor(data, output)` | B-factor profile PNG |
| `plot_ramachandran(data, output)` | Ramachandran scatter PNG |

### CLI
```bash
python3 summarize_structures.py 1abc.pdb
python3 summarize_structures.py Files/*.pdb --out summary.csv
python3 summarize_structures.py 1abc.pdb --missing
python3 summarize_structures.py 1abc.pdb --fasta [--chain A]
python3 summarize_structures.py 1abc.pdb --bfactor --plot
python3 summarize_structures.py 1abc.pdb --ramachandran --plot
python3 summarize_structures.py 1abc.pdb --rg
```

---

## Script: analyze_ligands.py

### Dependencies
BioPython, NumPy

### Key Functions
| Function | Returns | Description |
|---|---|---|
| `load_structure(path)` | `Structure` | Load PDB or mmCIF |
| `find_contacts(structure, sel1, sel2, cutoff)` | `list[Contact]` | Heavy-atom contacts |
| `find_hydrogen_bonds(...)` | `list[HBond]` | H-bonds via heavy-atom geometry |
| `find_pi_interactions(...)` | `list[PiInteraction]` | Pi-stacking + cation-pi |
| `find_salt_bridges(...)` | `list[SaltBridge]` | Charged group pairs |
| `analyze_ligand(structure, lig_sel, ...)` | `list[LigandReport]` | All interactions per ligand |

### CLI
```bash
python3 analyze_ligands.py structure.pdb --ligands
python3 analyze_ligands.py structure.pdb --contacts "resn ATP" "polymer" --cutoff 4.5
python3 analyze_ligands.py structure.pdb --hbonds "polymer" "resn ATP"
python3 analyze_ligands.py structure.pdb --pi "polymer" "resn ATP"
python3 analyze_ligands.py structure.pdb --salt-bridges "polymer" "polymer"
python3 analyze_ligands.py structure.pdb --analyze "resn ATP" --out report.csv
```

---

## Script: visualize_interactions.py

### Dependencies
BioPython, NumPy, Matplotlib; py3Dmol only needed for `build_view()` (notebook use)
HTML output loads **3Dmol.js via CDN** — no py3Dmol required for CLI use.

### Outputs
| File | Description |
|---|---|
| `{LIG}_{chain}{resi}.html` | Self-contained 3D viewer with per-type toggles and residue panel |
| `{LIG}_{chain}{resi}_fingerprint.png` | 2D dot-plot: residues × interaction types |

### HTML viewer features
- Toggle show/hide per interaction type (contacts, H-bonds, pi, salt bridges) — updates both 3D scene and residue panel simultaneously
- **Show all / Hide all** buttons
- Side panel lists every partner residue with atom names, distances, and roles
- Loads 3Dmol.js from CDN; fully self-contained single HTML file

### Color scheme
| Type | Color |
|---|---|
| Ligand | Yellow |
| Contacts | Cyan sticks |
| H-bonds | Orange sticks + yellow dashed lines |
| Pi interactions | Magenta sticks + magenta dashed lines |
| Salt bridges (+/−) | Blue/red sticks + orange dashed lines |

### CLI
```bash
python3 visualize_interactions.py structure.pdb --analyze "organic"
python3 visualize_interactions.py structure.pdb --analyze "resn ATP" --outdir ./images --width 1200
```

---

## Script: analyze_rmsd.py

### Dependencies
BioPython, NumPy, Matplotlib

### Key Functions
| Function | Description |
|---|---|
| `calculate_rmsd(path1, path2, ...)` | Cα RMSD after optimal superposition |
| `superpose(path1, path2, ..., output)` | Superpose mobile onto reference, save PDB |
| `per_residue_rmsd(path1, path2, ...)` | Per-residue Cα distance after superposition |
| `pairwise_rmsd_matrix(paths, ...)` | All-vs-all RMSD matrix |
| `nmr_ensemble_rmsd(path, ...)` | Per-model RMSD vs reference for NMR ensembles |

### CLI
```bash
python3 analyze_rmsd.py struct1.pdb struct2.pdb [--chain A] [--per-residue] [--plot]
python3 analyze_rmsd.py struct1.pdb struct2.pdb --superpose --out superposed.pdb
python3 analyze_rmsd.py Files/*.pdb --matrix --out rmsd.csv --plot
python3 analyze_rmsd.py nmr.pdb --ensemble --plot
```

---

## Script: ligand_rmsd.py

### Dependencies
BioPython, NumPy, tabulate

### Description
Computes heavy-atom RMSD for a named ligand across multiple PDB structures.
Aligns each structure to a reference by backbone Cα superposition (Kabsch).
Flags outliers by z-score and excludes them from the final average.

### Key Functions
| Function | Description |
|---|---|
| `get_backbone_atoms(structure, chain)` | Cα atoms for superposition |
| `get_ligand_heavy_atoms(structure, resn, chain)` | Non-hydrogen ligand atoms |
| `align_to_reference(mobile, ref, chain)` | In-place backbone superposition |
| `compute_centroid_coords(atom_dicts)` | Mean coords across all structures |
| `compute_rmsd(coords1, coords2)` | RMSD + atom count |

### CLI
```bash
python3 ligand_rmsd.py --ligand HEM Files/2HHB.pdb Files/4HHB.pdb
python3 ligand_rmsd.py --ligand HEM --reference ref.pdb "structures/*.pdb"
python3 ligand_rmsd.py --ligand HEM --chain A --zscore 1.5 "structures/*.pdb"
python3 ligand_rmsd.py --ligand HEM --no-exclude "structures/*.pdb"
```

### Test result (2HHB vs 4HHB, HEM)
- Centroid mode: both structures 0.126 Å from centroid (43 atoms matched)
- Explicit ref (4HHB): 2HHB deviates 0.251 Å — oxy vs deoxy heme rearrangement

---

## Script: conservation.py

### Dependencies
BioPython, NumPy, Matplotlib

### Description
Maps sequence conservation onto a PDB structure from a user-supplied multi-FASTA file.
No BLAST or network access required.

### Algorithm
1. Extract reference sequence from ATOM records
2. Pairwise-align each FASTA sequence to reference with BLOSUM62 global alignment
3. Compute per-column Shannon entropy → conservation score (1 = fully conserved, 0 = fully variable)
4. Write B-factor PDB, bar chart PNG, and CSV

### Outputs
| File | Description |
|---|---|
| `<stem>_conservation.png` | Bar chart, red→blue by conservation score |
| `<stem>_conservation.pdb` | PDB with conservation (0–100) in B-factor column |
| `<stem>_conservation.csv` | resnum, aa, score, entropy, gap_fraction, n_sequences |

### CLI
```bash
python3 conservation.py 1abc.pdb homologs.fasta
python3 conservation.py 1abc.pdb homologs.fasta --chain A --out results/
python3 conservation.py 1abc.pdb homologs.fasta --no-pdb
```

### PyMOL colouring
```
load 6D1Y_conservation.pdb;  spectrum b, red_white_blue, minimum=0, maximum=100
```

---

## Script: run_md.py

### Requirements
- **GROMACS ≥ 2019** (external, LGPL) — `brew install gromacs` / `sudo apt install gromacs`
- **numpy**, **matplotlib** — `pip install numpy matplotlib`
- **pdbfixer + openmm** (MIT) — strongly recommended; fills missing side-chain atoms
- **MDAnalysis** (GPL v2) — `pip install MDAnalysis` — enables DSSP and PCA

### Pipeline stages
| Stage | What it does |
|---|---|
| `prepare` | Strip HETATM/ANISOU; fix missing atoms via PDBFixer |
| `setup` | `pdb2gmx` → `editconf` (dodecahedral box) → `solvate` → `genion` |
| `minimize` | Steepest-descent EM (≤50,000 steps) |
| `equil` | NVT 100 ps (V-rescale thermostat) → NPT 100 ps (+ Berendsen barostat) |
| `run` | Production MD (Parrinello-Rahman barostat); XTC saved every 10 ps |
| `analyze` | RMSD, RMSF, Rg, energy, H-bonds, DSSP, PCA |
| `visualize` | matplotlib plots; PyMOL `.pml`; self-contained HTML report |

### Analysis outputs (in `<stem>_md/`)
| File | Description |
|---|---|
| `rmsd.xvg` | Backbone RMSD vs time |
| `rmsf.xvg` | Per-residue Cα RMSF |
| `gyrate.xvg` | Radius of gyration vs time |
| `energy.xvg` | Potential energy + temperature |
| `hbond.xvg` | Intra-protein H-bond count |
| `dssp.npz` | DSSP secondary structure array (frames × residues) |
| `pca.npz` | PCA results: variance, projected coordinates, PC1 eigenvector |
| `rmsf_bfactor.pdb` | Reference PDB with RMSF in B-factor column |
| `pc1_bfactor.pdb` | Reference PDB with PC1 eigenvector magnitude in B-factor column |
| `md_fit.xtc` | Trajectory centred and fitted to backbone |
| `plots/<stem>_summary.png` | All analysis panels in one figure |
| `plots/<stem>_rmsd.png` | Backbone RMSD time series |
| `plots/<stem>_rmsf.png` | Per-residue flexibility bar chart |
| `plots/<stem>_rg.png` | Radius of gyration time series |
| `plots/<stem>_energy.png` | Energy + temperature time series |
| `plots/<stem>_hbonds.png` | H-bond count time series |
| `plots/<stem>_dssp.png` | Secondary structure content stacked area |
| `plots/<stem>_pca.png` | PCA scree plot (variance per PC + cumulative) |
| `plots/<stem>_pca_scatter.png` | PC1 vs PC2 scatter coloured by time + time series |
| `plots/<stem>_pca_residues.png` | Per-residue PC1 eigenvector contribution |
| `<stem>_session.pml` | PyMOL session script (loads structure + trajectory) |
| `<stem>.pse` | PyMOL binary session |
| `<stem>_report.html` | Self-contained HTML report with embedded plots and statistics |

### Default simulation parameters
| Parameter | Value | Flag |
|---|---|---|
| Force field | AMBER99SB-ILDN | `--ff` |
| Water model | TIP3P | `--water` |
| Production length | 10 ns | `--ns` |
| Temperature | 300 K | `--temp` |
| Salt concentration | 0.15 M NaCl | `--conc` |
| Box padding | 1.0 nm | `--box` |
| Save interval | 10 ps | `--save-ps` |

### Known issues / fixes
- **Missing side chains**: Crystal structures with disordered residues cause `pdb2gmx` to fail ("atom CG not found"). Fix: install pdbfixer — it reconstructs missing atoms before `pdb2gmx` runs.
- **PyMOL trajectory**: `mset 1 xN` repeats state 1 N times (no movement). Use `mset 1 -N` to play states 1 through N sequentially. The script uses the correct form.
- **MDAnalysis TPR version**: MDAnalysis 2.x may not support the newest GROMACS TPR format. The pipeline falls back to `md_ref.pdb` (protein-only) as topology when TPR loading fails.
- **DSSP C-terminal atoms**: GROMACS names the C-terminal oxygens OC1/OC2. The pipeline renames them to O/OXT before running DSSP so backbone atom counts are equal.

### CLI
```bash
python3 run_md.py protein.pdb                            # full pipeline, 10 ns
python3 run_md.py protein.pdb --ns 100                   # 100 ns production
python3 run_md.py protein.pdb --steps prepare setup      # build inputs only
python3 run_md.py protein.pdb --steps analyze visualize  # post-process only
python3 run_md.py protein.pdb --ff charmm36m-iua --water tip4p
python3 run_md.py protein.pdb --ncores 8 --gpu
```

### Test run: 6D1Y at 10 ns (AMBER99SB-ILDN/TIP3P, 300 K)
- Backbone RMSD: 1.65 ± 0.26 Å (stable)
- PCA: PC1 = 35.8%, PC2 = 12.8%, cumulative = 48.5%
- DSSP: 101 frames × 320 residues
- All 20 output files generated successfully

---

## Git History
| Commit | Description |
|---|---|
| `125c467` | Initial commit: PDB downloader CLI |
| `790a576` | Add full PDB toolkit (analysis, MD pipeline) |
| `1025004` | Add ligand_rmsd.py |
| `846c83e` | Add ligand_rmsd examples to examples.md |
| `4b5c0f5` | Add requirements.txt |
| `021afe0` | Add interactive toggles and residue panel to HTML viewer |
