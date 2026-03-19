# PDB Toolkit — Project Plan

## Goal
A suite of standalone Python tools for downloading, analysing, visualising, and running molecular dynamics on protein structures from RCSB. Includes a Streamlit web app that wraps the three main analysis libraries.

## Project Structure
```
PDBs/
├── PLAN.md                    # This file — technical architecture
├── README.md                  # User-facing guide with examples
├── INSTALL.md                 # Setup instructions
├── examples.md                # Extended CLI examples
├── future_features.md         # Feature backlog
├── MD_TUTORIAL.md             # Science behind the MD pipeline
├── requirements.txt           # Python dependencies
├── app.py                     # Streamlit web app (5 tabs)
├── download_pdb.py            # Download PDB/CIF/FASTA from RCSB
├── summarize_structures.py    # Metadata, B-factor, Ramachandran, FASTA, Rg, BSA
├── analyze_ligands.py         # Ligand contacts, H-bonds, pi, salt bridges
├── visualize_interactions.py  # 3D HTML viewer + 2D LIGPLOT diagrams
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
streamlit run app.py          # web interface
```

---

## App: app.py

### Purpose
Streamlit browser-based interface wrapping `summarize_structures`, `analyze_ligands`, and `visualize_interactions`. No CLI knowledge required.

### Architecture
- **File persistence**: uploaded PDB saved to `tempfile.mkdtemp()` path stored in `st.session_state`; persists for the browser session
- **Analysis caching**: `@st.cache_data` keyed on `pdb_bytes: bytes` — identical file never re-analyzed; contact cutoff and BSA probe radius changes trigger re-analysis of their respective functions only
- **Structure objects**: `@st.cache_resource` keyed on path string for BioPython `Structure` (not picklable by `cache_data`)
- **Plots**: each plot function saves to a `NamedTemporaryFile`, bytes are read back and temp file deleted; bytes passed to `st.image()` and `st.download_button()`
- **3D viewer**: `_render_html()` called directly; resulting HTML string passed to `st.components.v1.html(height=720)`

### Tabs
| Tab | Key functions called |
|-----|---------------------|
| 📋 Structure Summary | `summarize_structure`, `bfactor_stats`, `plot_bfactor`, `ramachandran_data`, `plot_ramachandran`, `extract_fasta`, `find_missing_residues`, `get_ligands` |
| ⬛ Buried Surface Area | `buried_surface_areas`, `plot_bsa_matrix` |
| 🔬 Ligand Analysis | `primary_ligands`, `analyze_ligand` |
| 🌐 3D Viewer | `_render_html` (from `visualize_interactions`) |
| 🖼 2D Diagrams | `plot_interaction_summary`, `plot_ligand_2d` |

### Sidebar controls
- Contact cutoff (Å): 3.0–7.0, default 4.5 — triggers `_run_ligand_analysis` re-cache
- BSA probe radius (Å): 1.0–2.0, default 1.4 — triggers `_run_bsa` re-cache
- BSA sphere points: 50/100/200/500, default 100

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
BioPython ≥ 1.79, NumPy, Matplotlib

### Key Functions
| Function | Signature | Description |
|---|---|---|
| `load_structure` | `(path, sid) → Structure` | Load PDB or mmCIF |
| `parse_header` | `(path) → dict` | Resolution, R-work/free, method, organism, date |
| `summarize_structure` | `(path) → dict` | Full summary dict (all header + inventory fields) |
| `batch_summary` | `(paths, out_csv) → list[dict]` | Summarize N structures → CSV |
| `find_missing_residues` | `(path) → list[dict]` | SEQRES vs ATOM comparison per chain |
| `extract_fasta` | `(path, chain) → str` | FASTA from ATOM coordinates |
| `bfactor_stats` | `(structure, chain) → list[dict]` | Per-residue Cα B-factor with ±2σ flags |
| `radius_of_gyration` | `(structure, chain, atom_name) → float` | Cα Rg in Å |
| `ramachandran_data` | `(structure, chain) → list[dict]` | phi/psi angles + helix/sheet/allowed/outlier |
| `buried_surface_areas` | `(structure, probe_radius, n_points) → list[dict]` | Pairwise chain BSA via Shrake-Rupley; sorted by bsa_total desc |
| `plot_bfactor` | `(data, output, title) → Path` | B-factor profile PNG |
| `plot_ramachandran` | `(data, output, title) → Path` | Ramachandran scatter PNG |
| `plot_bsa_matrix` | `(bsa_results, chain_ids, output, title) → Path` | BSA heatmap + ranked bar chart PNG |

### BSA algorithm
```
BSA(A, B) = ( SASA_A + SASA_B − SASA_AB ) / 2
bsa_on_A  = SASA_A − SASA_A_in_complex
bsa_on_B  = SASA_B − SASA_B_in_complex
```
- Uses `Bio.PDB.SASA.ShrakeRupley(probe_radius, n_points)`
- Per-chain SASA within the complex obtained by summing per-atom SASA for atoms belonging to each chain after running on the combined sub-structure
- Helper `_chain_substructure(source, chain_ids)` builds a new Structure from copied chains

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

Note: `buried_surface_areas` is importable but not yet exposed as a CLI flag. Use `app.py` or import directly.

---

## Script: analyze_ligands.py

### Dependencies
BioPython, NumPy

### Key Functions
| Function | Returns | Description |
|---|---|---|
| `load_structure(path)` | `Structure` | Load PDB or mmCIF; fills missing element fields |
| `get_ligands(structure)` | `list[Residue]` | All non-water HETATM residues |
| `primary_ligands(structure, min_atoms, exclude_excipients)` | `list[Residue]` | Non-excipient ligands ≥ min_atoms heavy atoms, sorted largest first |
| `find_contacts(structure, sel1, sel2, cutoff)` | `list[Contact]` | Heavy-atom contacts |
| `find_hydrogen_bonds(structure, sel1, sel2, ...)` | `list[HBond]` | H-bonds via heavy-atom geometry (N/O/S); angle filter |
| `find_pi_interactions(structure, sel1, sel2, pdb_path, ...)` | `list[PiInteraction]` | Pi-stacking + cation-pi; ligand rings from CONECT records |
| `find_salt_bridges(structure, sel1, sel2, cutoff)` | `list[SaltBridge]` | ARG/LYS/HIS ↔ ASP/GLU charged pairs |
| `analyze_ligand(structure, lig_sel, protein_sel, ...)` | `list[LigandReport]` | All interactions per ligand; runs all four above |

### Constants
- `WATER_NAMES` — frozenset of water residue names (HOH, WAT, H2O, DOD, SOL)
- `COMMON_EXCIPIENTS` — frozenset of ~50 crystallographic additives filtered by `primary_ligands()`
- `PROTEIN_RESIDUES`, `NUCLEIC_RESIDUES` — frozensets for polymer classification
- `PROTEIN_AROMATIC_RINGS` — dict mapping PHE/TYR/TRP/HIS to ring atom name lists
- `CATION_ATOMS`, `ANION_ATOMS` — dict mapping residue names to charged atom names
- Default cutoffs: contacts 4.5 Å, H-bonds 3.5 Å, salt bridges 4.0 Å, pi 5.5 Å, cation-pi 6.0 Å

### Data classes
- `Contact` — chain1, resn1, resi1, atom1, chain2, resn2, resi2, atom2, distance
- `HBond` — donor_chain/resn/resi/atom, acceptor_chain/resn/resi/atom, distance, angle
- `PiInteraction` — chain1/resn1/resi1/ring1_label, chain2/resn2/resi2/ring2_label, center_distance, plane_angle, subtype (face_to_face/edge_to_face/intermediate/cation_pi)
- `SaltBridge` — cation_chain/resn/resi/atom, anion_chain/resn/resi/atom, distance
- `LigandReport` — aggregates all four for one ligand; properties: n_contacts, n_hbonds, n_pi, n_salt_bridges

### Selection parser
`SelectionParser(structure)` — PyMOL-compatible subset:
- Keywords: `all`, `none`, `polymer` (`.protein`/`.nucleic`), `hetatm`, `organic`, `solvent`, `water`, `not`, `and`, `or`, `within N of`
- Properties: `chain`, `resn`, `resi` (ranges, `+` lists), `name`, `elem`, `b` (comparisons)
- Digit-prefixed residue names (e.g. `1PE`, `2HB`) handled by token merging

### CLI
```bash
python3 analyze_ligands.py structure.pdb              # auto-detect + full report
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
BioPython, NumPy, Matplotlib; `py3Dmol` only for `build_view()` (notebook use).
HTML output loads **3Dmol.js via CDN** — py3Dmol not required for CLI use.

### Key Functions
| Function | Description |
|---|---|
| `save_html(structure_path, report, output_path, ...)` | Write self-contained 3D HTML viewer |
| `_render_html(structure_path, report, structure, width, height, active_types)` | Generate HTML string (called by `save_html`; used directly in `app.py`) |
| `plot_interaction_summary(report, output_path, types)` | 2D fingerprint dot-plot PNG |
| `plot_ligand_2d(structure_path, report, output_path, structure)` | LIGPLOT-style 2D diagram PNG |
| `save_type_images(structure_path, report, outdir, stem, ...)` | Per-type HTML + filtered fingerprint for each type with data |
| `visualize(structure_path, ligand_sel, protein_sel, outdir, ...)` | High-level: analysis + all outputs |

### Outputs per ligand
| File | Description |
|---|---|
| `{LIG}_{chain}{resi}.html` | Self-contained 3D viewer with per-type toggles and residue panel |
| `{LIG}_{chain}{resi}_fingerprint.png` | 2D dot-plot: residues × interaction types |
| `{LIG}_{chain}{resi}_2d.png` | LIGPLOT-style 2D diagram |
| `{LIG}_{chain}{resi}_{type}.html` | Per-type 3D viewer (only that interaction pre-checked) |
| `{LIG}_{chain}{resi}_{type}_fingerprint.png` | Per-type filtered fingerprint |

### HTML viewer internals
- PDB content embedded as a JS string literal; no server needed
- `updateScene()` JS function called on checkbox change: resets base cartoon, re-applies atom styles in priority order (contacts < pi < salt < hbonds), removes/adds shapes and labels
- All atom styles, dashed-line cylinder specs, and label specs pre-serialised to JSON by `_render_html()` and injected into the JS `atomStyles`, `shapeSpecs`, `labelSpecs` variables

### CLI
```bash
python3 visualize_interactions.py structure.pdb                       # auto-detect
python3 visualize_interactions.py structure.pdb --analyze "resn ATP"
python3 visualize_interactions.py structure.pdb --analyze "resn ATP" --outdir ./images --width 1200
```

---

## Script: analyze_rmsd.py

### Dependencies
BioPython, NumPy, Matplotlib

### Key Functions
| Function | Description |
|---|---|
| `calculate_rmsd(path1, path2, chain, align)` | Cα RMSD after optimal superposition |
| `superpose(path1, path2, chain, align, output)` | Superpose mobile onto reference, save PDB |
| `per_residue_rmsd(path1, path2, chain, align, threshold)` | Per-residue Cα distance after superposition |
| `pairwise_rmsd_matrix(paths, chain, align)` | All-vs-all RMSD matrix |
| `nmr_ensemble_rmsd(path, chain, ref_model)` | Per-model RMSD vs reference for NMR ensembles |

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
- **GROMACS ≥ 2019** — `brew install gromacs` / `sudo apt install gromacs`
- **numpy**, **matplotlib**
- **pdbfixer + openmm** — fills missing side-chain atoms (strongly recommended)
- **MDAnalysis** — enables DSSP timeline and PCA

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
- **Missing side chains**: Use pdbfixer — it reconstructs missing atoms before `pdb2gmx` runs.
- **PyMOL trajectory**: `mset 1 xN` (wrong) repeats state 1; use `mset 1 -N` to play states 1–N. Script uses the correct form.
- **MDAnalysis TPR version**: Falls back to `md_ref.pdb` as topology when TPR loading fails.
- **DSSP C-terminal atoms**: Pipeline renames GROMACS OC1/OC2 to O/OXT before DSSP so backbone atom counts match.

### CLI
```bash
python3 run_md.py protein.pdb                            # full pipeline, 10 ns
python3 run_md.py protein.pdb --ns 100
python3 run_md.py protein.pdb --steps prepare setup
python3 run_md.py protein.pdb --steps analyze visualize
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
| `9ac6b1d` | Add LIGPLOT-style 2D diagram, per-type images, smart ligand auto-detection |
| `6a090a2` | Update README.md |
| `c04630f` | Fix element guessing, digit-prefixed resn parsing, CSV pdb_id truncation |
| *(current)* | Add Streamlit app (app.py), buried surface area (buried_surface_areas, plot_bsa_matrix), update all markdown files |
