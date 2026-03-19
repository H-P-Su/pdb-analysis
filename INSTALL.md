# Installation

## Requirements

- Python 3.9 or later
- Git (to clone the repository)
- GROMACS ≥ 2019 — only required for `run_md.py`

---

## Quick setup

```bash
git clone https://github.com/YOUR_USERNAME/pdb-toolkit.git
cd pdb-toolkit

python3 -m venv venv
source venv/bin/activate          # macOS / Linux
# venv\Scripts\activate.bat       # Windows

pip install -r requirements.txt
```

### Verify

```bash
python3 download_pdb.py 1TIM
# → Files/1TIM.pdb should appear

python3 summarize_structures.py Files/1TIM.pdb
# → prints metadata table

streamlit run app.py
# → opens browser at http://localhost:8501
```

---

## Core dependencies (`requirements.txt`)

| Package | Version | Used by |
|---------|---------|---------|
| `biopython` | ≥ 1.79 | all analysis scripts |
| `numpy` | ≥ 1.24 | all analysis scripts |
| `matplotlib` | ≥ 3.7 | plots in summarize, visualize, rmsd, conservation |
| `scipy` | ≥ 1.10 | statistical utilities |
| `tabulate` | ≥ 0.9 | `ligand_rmsd.py` table output |
| `py3Dmol` | ≥ 2.0 | `build_view()` in notebooks (HTML output does not require it) |
| `streamlit` | ≥ 1.30 | `app.py` web interface |
| `pandas` | ≥ 2.0 | `app.py` data tables |

Install all at once:

```bash
pip install -r requirements.txt
```

---

## Optional: Streamlit web app

`app.py` is a browser-based interface that wraps `summarize_structures`, `analyze_ligands`, and `visualize_interactions` into four tabs — no command line needed.

```bash
streamlit run app.py
# Open http://localhost:8501 in any browser
```

Upload any `.pdb` file via the sidebar. All analyses run automatically and results are downloadable as CSV or PNG.

To run on a specific port or expose to your network:

```bash
streamlit run app.py --server.port 8888
streamlit run app.py --server.address 0.0.0.0   # accessible from other machines
```

---

## Optional: Molecular dynamics pipeline

`run_md.py` requires GROMACS as an external binary:

```bash
brew install gromacs              # macOS (Homebrew)
sudo apt install gromacs          # Ubuntu / Debian
conda install -c conda-forge gromacs   # Conda
module load gromacs               # HPC cluster — check your module name
```

Two Python packages are also required inside the venv:

```bash
pip install pdbfixer openmm       # fills missing side-chain atoms
pip install MDAnalysis            # DSSP secondary structure + PCA
```

**Why pdbfixer?** Most crystal structures have disordered loop residues with missing side-chain heavy atoms. GROMACS `pdb2gmx` fails if any heavy atom is absent. PDBFixer reconstructs them before topology generation.

**Why MDAnalysis?** Enables the DSSP per-frame secondary structure timeline and backbone PCA in the `analyze` stage. Without it those two analyses are silently skipped and the summary plot has fewer panels.

### Verify GROMACS

```bash
gmx --version
# Should print GROMACS version ≥ 2019
```

---

## Platform notes

| Platform | Notes |
|----------|-------|
| macOS (Apple Silicon) | All packages including BioPython install natively via pip. GROMACS from Homebrew works. |
| macOS (Intel) | Same as above. |
| Ubuntu / Debian | Standard `apt` GROMACS build is usually sufficient. GPU builds require the CUDA-enabled package from GROMACS website. |
| Windows | Use WSL2 (Ubuntu) for the full toolkit. Native Windows Python works for all non-MD scripts. |
| HPC cluster | Load a GROMACS module, then install Python packages into a user venv with `pip install --user`. |

---

## File layout after setup

```
PDBs/
├── venv/                  # virtual environment (not in git)
├── Files/                 # downloaded structures (not in git)
├── demo_output/           # example outputs (tracked for reference)
├── app.py                 # Streamlit web interface
├── download_pdb.py
├── summarize_structures.py
├── analyze_ligands.py
├── visualize_interactions.py
├── analyze_rmsd.py
├── ligand_rmsd.py
├── conservation.py
├── run_md.py
├── examples.py
├── requirements.txt
├── README.md
├── PLAN.md
├── INSTALL.md             # this file
├── examples.md
├── future_features.md
└── MD_TUTORIAL.md
```
