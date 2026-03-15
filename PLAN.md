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

## Dependencies

Install into the venv:
```bash
python3 -m venv venv && source venv/bin/activate
pip install biopython numpy
```

---

## Script: analyze_ligands.py

### Requirements
- Depends on **BioPython** and **NumPy** (install into venv)
- All public functions importable for use in other projects
- PyMOL-compatible selection syntax (subset) for specifying atom groups
- Interaction types: contacts, H-bonds, pi interactions, salt bridges
- Ligand aromatic ring detection from CONECT records (PDB only)

### Importable API
| Function | Returns | Description |
|---|---|---|
| `load_structure(path)` | `Structure` | Load PDB or mmCIF |
| `get_ligands(structure)` | `list[Residue]` | All non-water HETATM residues |
| `find_contacts(structure, sel1, sel2, cutoff)` | `list[Contact]` | Heavy-atom contacts |
| `find_hydrogen_bonds(structure, sel1, sel2, ...)` | `list[HBond]` | H-bonds via heavy-atom geometry |
| `find_pi_interactions(structure, sel1, sel2, ...)` | `list[PiInteraction]` | Pi-stacking + cation-pi |
| `find_salt_bridges(structure, sel1, sel2, ...)` | `list[SaltBridge]` | Charged group pairs |
| `analyze_ligand(structure, lig_sel, ...)` | `list[LigandReport]` | All interactions per ligand |
| `to_csv(interactions, path)` | — | Write results to CSV |
| `print_table(interactions)` | — | Pretty-print to terminal |

### Selection syntax (PyMOL-compatible)
- `all`, `none`, `polymer`, `hetatm`, `organic`, `solvent`
- `chain A+B`, `resn LIG+ATP`, `resi 100`, `resi 100-200`, `name CA`
- `elem N+O`, `b > 50`
- `not`, `and`, `or`, `( )`
- `within N of SELECTION`

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

---

## Script: visualize_interactions.py

### Requirements
- Depends on `analyze_ligands.py`, **py3Dmol**, **matplotlib**, **numpy**
- Generates two outputs per ligand: interactive HTML + 2D fingerprint PNG
- All functions importable

### Outputs
| File | Description |
|---|---|
| `{LIG}_{chain}{resi}.html` | Interactive 3D viewer (py3Dmol, opens in any browser) |
| `{LIG}_{chain}{resi}_fingerprint.png` | 2D dot-plot: residues × interaction types |

### Color scheme
| Type | Color |
|---|---|
| Ligand | Yellow |
| Contacts | Cyan sticks |
| H-bonds | Orange sticks + yellow dashed lines |
| Pi interactions | Magenta sticks + magenta dashed lines |
| Salt bridges (+/−) | Blue/red sticks + orange dashed lines |

### Importable API
| Function | Description |
|---|---|
| `build_view(path, report, structure, w, h)` | Build py3Dmol view (for notebooks) |
| `save_html(path, report, output, ...)` | Write interactive HTML file |
| `plot_interaction_summary(report, output)` | Write 2D fingerprint PNG |
| `visualize(path, lig_sel, ...)` | Run full pipeline, return list of outputs |

### CLI
```bash
python3 visualize_interactions.py structure.pdb --analyze "organic"
python3 visualize_interactions.py structure.pdb --analyze "resn ATP" --outdir ./images --width 1200
```

---

## Script: summarize_structures.py

### Requirements
- Depends on **BioPython**, **NumPy**, **Matplotlib**
- All functions importable

### Functions
| Function | Description |
|---|---|
| `parse_header(path)` | Resolution, R-work/free, method, organism, date from REMARK records |
| `summarize_structure(path)` | Full summary dict (chains, residues, atoms, ligands, missing) |
| `batch_summary(paths, out_csv)` | Summarize N structures → CSV |
| `find_missing_residues(path)` | SEQRES vs ATOM comparison per chain |
| `extract_fasta(path, chain)` | FASTA from ATOM coordinates |
| `bfactor_stats(structure, chain)` | Per-residue Cα B-factor with ±2σ flags |
| `radius_of_gyration(structure, chain)` | Cα Rg in Å |
| `ramachandran_data(structure, chain)` | phi/psi angles + region classification |
| `plot_bfactor(data, output)` | B-factor profile bar plot PNG |
| `plot_ramachandran(data, output)` | Ramachandran scatter plot PNG |

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

## Script: analyze_rmsd.py

### Requirements
- Depends on **BioPython** (`Superimposer`, `PairwiseAligner`), **NumPy**, **Matplotlib**
- All functions importable

### Functions
| Function | Description |
|---|---|
| `load_ca_atoms(path, chain, model)` | Cα atoms + one-letter sequence |
| `align_sequences(seq1, seq2)` | Global alignment → matched index pairs |
| `calculate_rmsd(path1, path2, ...)` | Cα RMSD after optimal superposition |
| `superpose(path1, path2, ..., output)` | Superpose mobile onto reference, save PDB |
| `per_residue_rmsd(path1, path2, ...)` | Per-residue Cα distance after superposition |
| `pairwise_rmsd_matrix(paths, ...)` | All-vs-all RMSD matrix |
| `nmr_ensemble_rmsd(path, ...)` | Per-model RMSD vs reference for NMR ensembles |
| `plot_rmsd_matrix(matrix, labels, out)` | Heatmap PNG |
| `plot_per_residue_rmsd(data, out)` | Bar plot PNG |
| `plot_ensemble_rmsd(data, out)` | Ensemble bar plot PNG |

### Alignment modes
- `align_by="sequence"` — global pairwise alignment (handles insertions/deletions)
- `align_by="resnum"` — match by residue number (fast; requires same numbering scheme)

### CLI
```bash
python3 analyze_rmsd.py struct1.pdb struct2.pdb
python3 analyze_rmsd.py struct1.pdb struct2.pdb --chain A --per-residue --plot
python3 analyze_rmsd.py struct1.pdb struct2.pdb --superpose --out superposed.pdb
python3 analyze_rmsd.py Files/*.pdb --matrix --out rmsd.csv --plot
python3 analyze_rmsd.py nmr.pdb --ensemble --plot
```

---

## Script: run_md.py

### Requirements
- Requires **GROMACS ≥ 2019** (external executable, LGPL v2.1) — `brew install gromacs` / `sudo apt install gromacs`
- **numpy**, **matplotlib** (BSD) — `pip install numpy matplotlib`
- **pdbfixer + openmm** (MIT) — `pip install pdbfixer openmm` — **strongly recommended**; fills missing side-chain atoms that would otherwise cause `pdb2gmx` to fail (e.g. disordered residues in crystal structures)
- **MDAnalysis** (GPL v2) — `pip install MDAnalysis` — recommended; enables DSSP secondary structure analysis

### Pipeline stages
| Stage | What it does |
|---|---|
| `prepare` | Strip HETATM/ANISOU; fix missing atoms via PDBFixer if available |
| `setup` | `pdb2gmx` (force field + H atoms) → `editconf` (dodecahedral box) → `solvate` (TIP3P water) → `genion` (Na⁺/Cl⁻ neutralisation) |
| `minimize` | Steepest-descent energy minimisation (≤50,000 steps, F_max < 1000 kJ/mol/nm) |
| `equil` | NVT 100 ps (V-rescale thermostat, position restraints) → NPT 100 ps (+ Berendsen barostat) |
| `run` | Production MD with Parrinello-Rahman barostat; XTC saved every 10 ps |
| `analyze` | Backbone RMSD, per-residue Cα RMSF, radius of gyration, potential energy + temperature, H-bond count, DSSP secondary structure timeline |
| `visualize` | matplotlib multi-panel summary PNG + individual plots; PyMOL `.pml` session script with trajectory loaded as states, RMSF B-factor colouring, and `mset 1 -N` movie setup |

### Output files (in `<stem>_md/`)
| File | Description |
|---|---|
| `topol.top` | GROMACS topology (force-field parameters) |
| `em.gro` | Energy-minimised coordinates |
| `md.tpr` | Production run binary (topology + coordinates) |
| `md.xtc` | Raw trajectory (all atoms) |
| `md_fit.xtc` | Trajectory centred and fitted to backbone |
| `rmsd.xvg` | Backbone RMSD vs time |
| `rmsf.xvg` | Per-residue Cα RMSF |
| `gyrate.xvg` | Radius of gyration vs time |
| `energy.xvg` | Potential energy + temperature vs time |
| `hbond.xvg` | Intra-protein H-bond count vs time |
| `dssp.npz` | DSSP secondary structure array (frames × residues) |
| `rmsf_bfactor.pdb` | Reference PDB with RMSF in B-factor column |
| `plots/<stem>_summary.png` | All analysis panels in one figure |
| `<stem>_session.pml` | PyMOL session script |
| `<stem>.pse` | PyMOL binary session (saved by the .pml) |

### Default simulation parameters
| Parameter | Value | Flag to change |
|---|---|---|
| Force field | AMBER99SB-ILDN | `--ff` |
| Water model | TIP3P | `--water` |
| Production length | 10 ns | `--ns` |
| Temperature | 300 K | `--temp` |
| Salt concentration | 0.15 M NaCl | `--conc` |
| Box padding | 1.0 nm | `--box` |
| Save interval | 10 ps | `--save-ps` |

### Known issues / fixes
- Crystal structures with disordered side chains (common) will cause `pdb2gmx` to fail with "atom CG not found". Fix: install pdbfixer — it reconstructs missing atoms before `pdb2gmx` runs.
- PyMOL trajectory playback: `mset 1 xN` repeats state 1 N times (no movement). The correct syntax is `mset 1 -N` which plays states 1 through N sequentially. The script uses the correct form.

### CLI
```bash
python3 run_md.py protein.pdb                            # full pipeline, 10 ns
python3 run_md.py protein.pdb --ns 100                   # 100 ns production
python3 run_md.py protein.pdb --steps prepare setup      # build inputs only
python3 run_md.py protein.pdb --steps analyze visualize  # post-process only
python3 run_md.py protein.pdb --ff charmm36m-iua --water tip4p
python3 run_md.py protein.pdb --ncores 8 --gpu           # parallel/GPU run
python3 run_md.py protein.pdb --ns 0.1 --nvt-ps 10 --npt-ps 10  # quick test run
```

---

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
