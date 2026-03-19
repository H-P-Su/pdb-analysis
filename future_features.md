# Future Features

Ideas for extending this project. Checked items are fully implemented.

---

## Structure Inventory / Parsing

- [x] Parse and summarize a PDB file: chain count, residue count, atom count, ligands present
- [x] Detect missing residues or sequence gaps (compare SEQRES vs ATOM records)
- [x] Extract FASTA sequence directly from ATOM coordinates (independent of download)
- [x] List all ligands/small molecules (HETATM records) with residue names and chain
- [ ] Extract individual chains to separate PDB files
- [ ] Strip solvent/ligands — output protein-only or nucleic-acid-only PDB
- [ ] Extract a specific HETATM ligand to its own file (for docking prep)
- [ ] PDB → mmCIF conversion locally (without re-downloading)

---

## Geometry & Conformation

- [x] **Ramachandran plot** — phi/psi backbone dihedrals; outlier flagging; PNG output
- [x] **B-factor analysis** — per-residue Cα B-factor; ±2σ flagging; profile plot PNG
- [x] **Radius of gyration** — Cα Rg in Å (or all heavy atoms)
- [x] **Secondary structure assignment** — DSSP per-residue per-frame via MDAnalysis; stacked area chart over MD trajectory
- [x] **Buried surface area** — Shrake-Rupley SASA for every chain pair; BSA heatmap matrix; per-chain interface %; `buried_surface_areas()` + `plot_bsa_matrix()` in `summarize_structures.py`; ⬛ BSA tab in `app.py`
- [ ] **Half-sphere exposure (HSE)** — `Bio.PDB.HSExposure.HSExposureCB`; better burial metric than raw SASA for identifying packing; complement to B-factor
- [ ] **Residue depth** — `Bio.PDB.ResidueDepth`; distance of each residue Cα from the molecular surface; useful for identifying buried vs. exposed residues
- [ ] **Disulfide bond detection** — find CYS–CYS pairs within ~2.05 Å S–S distance

---

## Interactions

- [x] **Salt bridges** — protein–protein and protein–ligand (ARG/LYS/HIS vs ASP/GLU)
- [x] **Hydrogen bond analysis** — heavy-atom geometry with angle filter; no explicit H needed
- [x] **Pi interactions** — face-to-face, edge-to-face stacking + cation-pi (ARG/LYS/HIS)
- [x] **Ligand contacts** — all residues within N Å of each ligand (default 4.5 Å)
- [ ] **Disulfide bonds** — CYS SG–SG distance ≤ 2.1 Å
- [ ] **Hydrophobic core identification** — cluster buried nonpolar residues (LEU/ILE/VAL/PHE/MET/ALA); use residue depth or HSE for burial filter
- [ ] **Steric clash detection** — flag heavy-atom pairs closer than 0.6 × (sum of van der Waals radii); `NeighborSearch` already in place

---

## Sequence & Conservation

- [x] **Sequence conservation mapping** — per-residue Shannon entropy from multi-FASTA; BLOSUM62 pairwise alignment; outputs bar chart PNG, B-factor PDB, CSV
- [ ] **Sequence entropy heatmap** — multi-chain conservation side by side in one figure

---

## Multi-structure / Comparative

- [x] **RMSD** — Cα RMSD between two structures; auto sequence alignment; resnum fallback
- [x] **Structure superposition** — superpose mobile onto reference, save transformed PDB
- [x] **Pairwise RMSD matrix** — all-vs-all for N structures; heatmap PNG + CSV export
- [x] **Per-residue RMSD** — per-residue Cα distance after superposition; bar chart PNG
- [x] **NMR ensemble analysis** — per-model RMSD vs reference model; bar chart PNG
- [ ] **CE alignment** — `Bio.PDB.cealign.CEAligner`; sequence-independent structural alignment; useful for distantly related structures where sequence alignment fails
- [ ] **Structural diff** — highlight residues that differ significantly between two structures (RMSD > threshold); output annotated PDB

---

## Structural Validation / QC

- [x] **Header metadata extraction** — resolution, R-factor, deposition date, organism
- [x] **Completeness report** — SEQRES vs ATOM residue count, % modeled per chain
- [ ] **Steric clash detection** — flag atom pairs below van der Waals contact threshold
- [ ] **Bond length / angle outliers** — flag covalent geometry deviations
- [ ] **Rotamer outlier detection** — compare sidechain χ angles to library (MolProbity-style)

---

## Visualization

- [x] **Interactive 3D HTML viewer** — 3Dmol.js (CDN); per-type toggle checkboxes (contacts/H-bonds/pi/salt); Show all / Hide all; side panel with atom names and distances; per-type separate HTMLs
- [x] **2D interaction fingerprint** — dot-plot (residues × interaction types) PNG
- [x] **LIGPLOT-style 2D diagram** — SVD-projected ligand; CPK atoms; spoked arcs; no RDKit needed
- [x] **BSA heatmap matrix** — chains × chains heatmap + ranked bar chart
- [x] **PyMOL .pml script generation** — session script loading structure + XTC trajectory; RMSF + PC1 B-factor colouring
- [ ] **Binding site surface mesh** — render molecular surface for the pocket using `Bio.PDB.get_surface()` or MSMS
- [ ] **Conservation-coloured 3D viewer** — embed conservation scores into the HTML viewer as a colour dimension

---

## Web Interface (`app.py`)

- [x] **Streamlit app** — file upload → full analysis; 5 tabs: Structure Summary, BSA, Ligand Analysis, 3D Viewer, 2D Diagrams
- [x] **Embedded 3Dmol.js viewer** — rendered inside Streamlit via `st.components.v1.html`
- [x] **Downloadable outputs** — every PNG and CSV downloadable from the UI
- [ ] **Multi-file upload** — upload several PDBs, get batch summary table + pairwise RMSD matrix
- [ ] **Conservation tab** — upload FASTA alongside PDB, render conservation bar chart in app
- [ ] **RMSD comparison tab** — upload two PDBs, show per-residue RMSD bar chart + superposed 3D viewer

---

## Molecular Dynamics

- [x] **Full MD pipeline** — GROMACS: prep → topology → solvate → minimise → NVT/NPT equil → production → analysis → plots + PyMOL session (`run_md.py`)
- [x] **MD trajectory analysis** — RMSD, per-residue RMSF, radius of gyration, energy/temperature, H-bond count, DSSP
- [x] **Principal component analysis** — backbone PCA; scree plot, PC1 vs PC2 scatter, per-residue PC1 eigenvector; PC1 B-factor PDB
- [x] **Automated HTML report** — self-contained single-file report with embedded plots, statistics, parameter table, PyMOL instructions
- [x] **Ligand RMSD across structures** — backbone Cα superposition + ligand heavy-atom RMSD; centroid or explicit reference; z-score outlier exclusion
- [ ] **Ligand MD prep** — ACPYPE / ATB integration for small-molecule force-field parameters; merge ligand topology with protein topology
- [ ] **Free energy perturbation** — alchemical FEP setup for relative binding affinity calculations
- [ ] **Binding site volume tracking** — monitor pocket volume over MD trajectory (fpocket / MDAnalysis `Hole`)
- [ ] **Protein–protein interface tracking** — BSA vs time over MD trajectory

---

## Batch / Workflow

- [x] **Batch structure summary** — CSV with resolution, R-factors, chains, residues, ligands, organism, missing residues
- [ ] **Watch a list** — poll RCSB for updates to a set of IDs and re-download if revised
- [ ] **Pipeline YAML** — define a full analysis workflow as YAML; run with a single command
