# Molecular Dynamics Tutorial

This tutorial explains what the `run_md.py` pipeline does at each stage — the science behind the choices, what each GROMACS tool produces, and how to interpret the outputs. A practical quick-start guide for running a new structure is at the end.

---

## What is Molecular Dynamics?

A molecular dynamics (MD) simulation numerically integrates Newton's equations of motion for every atom in the system. Given atom positions and a mathematical description of interatomic forces (the **force field**), the simulation advances time in discrete steps (typically 2 femtoseconds), computing new velocities and positions at each step.

Over millions of steps, this produces a **trajectory** — a movie of atomic motion at physiological temperature and pressure. Unlike a static crystal structure, an MD trajectory reveals:

- Which regions of the protein are flexible vs. rigid
- How secondary structure is maintained or lost over time
- How a protein breathes, opens, closes, or unfolds
- Whether a ligand-bound structure is stable in solution
- The distribution of conformations a protein samples

The pipeline in `run_md.py` automates the seven stages needed to go from a PDB file to an analysed, visualised trajectory:

```
PDB file
   │
   ▼  prepare    Fix missing atoms, strip ligands/water
   │
   ▼  setup      Build force-field topology, add water box and ions
   │
   ▼  minimize   Remove steric clashes (energy minimisation)
   │
   ▼  equil      Gradually heat and equilibrate (NVT → NPT)
   │
   ▼  run        Production MD — the actual simulation
   │
   ▼  analyze    RMSD, RMSF, Rg, energy, H-bonds, secondary structure
   │
   ▼  visualize  Time-series plots (PNG) + PyMOL session (.pml)
```

---

## Stage 1 — Prepare (`prepare_structure`)

### What the code does

```python
fixer = PDBFixer(str(cfg.pdb_in))
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
```

### Why this is necessary

Crystal structures deposited in the PDB are derived from electron density maps. When a side chain is disordered — pointing into solvent and adopting multiple conformations — crystallographers often omit those atoms rather than model them inaccurately. A residue like LYS might be present in the SEQRES record but missing its CG, CD, CE, or NZ atoms in the ATOM records.

`pdb2gmx` (the GROMACS topology builder) requires every atom listed in its force-field template for a residue to be present. If even one side-chain atom is missing, it exits with a fatal error:

```
Fatal error: atom CG in residue LYS 2 not found
```

PDBFixer reconstructs missing atoms by placing them at geometrically reasonable positions based on ideal bond lengths and angles, using rotamer libraries for side chains. This gives `pdb2gmx` a complete structure to work from.

**Why ligands are removed:** Standard force fields (AMBER, CHARMM, OPLS) contain parameters only for the 20 canonical amino acids, DNA/RNA nucleotides, and common ions. Small-molecule ligands need their own parameter sets, generated separately with tools like ACPYPE or the ATB server. Mixing an unknown ligand into the topology without parameters causes `pdb2gmx` to fail. The pipeline removes all HETATM records by default and simulates the apo protein.

**Fallback without PDBFixer:** A basic line filter strips HETATM and ANISOU records. This works only if the protein has no missing heavy atoms — uncommon for real crystal structures.

---

## Stage 2 — Setup (`setup_topology`)

Four sub-steps build the full simulation system from the cleaned PDB.

### 2a. `pdb2gmx` — force-field assignment

```bash
gmx pdb2gmx -f prepared.pdb -o processed.gro -p topol.top \
            -ff amber99sb-ildn -water tip3p -ignh
```

**What it does:**
- Reads the sequence and assigns AMBER99SB-ILDN atom types to every atom
- Deletes all existing hydrogen atoms (`-ignh`) and adds them back according to the force-field rules (correct protonation at pH 7)
- Writes `topol.top` — the topology file describing every bond, angle, dihedral, and non-bonded interaction parameter in the system
- Writes `posre.itp` — position restraint parameters for protein heavy atoms, used during equilibration

**Why AMBER99SB-ILDN?** This force field was developed specifically for protein backbone and side-chain dynamics. The "ILDN" suffix refers to improved parameters for isoleucine (I), leucine (L), aspartate (D), and asparagine (N) rotamers, validated against NMR order parameters. It is the most widely benchmarked force field for globular proteins.

**Output file format (GRO):** GROMACS uses its own coordinate format (`.gro`) internally. Each line encodes residue name, atom name, atom number, x/y/z position in nanometres, and optionally velocity.

### 2b. `editconf` — periodic box

```bash
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt dodecahedron
```

MD simulations use **periodic boundary conditions (PBC)**: the simulation box is replicated in all directions, so atoms exiting one face immediately re-enter from the opposite face. This eliminates surface effects and allows a small box to represent bulk solution.

**Box shape — dodecahedron:** A rhombic dodecahedron is the most efficient shape for a roughly spherical protein. It has ~71% the volume of a cubic box of the same minimum image distance, meaning ~29% fewer water molecules to simulate — a significant speed-up.

**`-d 1.0`:** The protein surface is at least 1.0 nm from the box edge in every direction. This ensures that the protein never interacts with its own periodic image across the box boundary (the minimum image convention requires the interaction cutoff — 1.0 nm — to be less than half the box dimension).

### 2c. `solvate` — add water

```bash
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
```

The empty space in the box is filled with pre-equilibrated TIP3P water molecules taken from a reference configuration (`spc216.gro`). Water molecules that overlap with the protein are removed. The topology is updated with the correct number of water molecules.

**Why explicit water?** Implicit solvent models approximate water as a continuum, which is fast but misses hydrogen bonding, hydrophobic effects driven by water structure, and the kinetics of water exchange around the protein. Explicit TIP3P water, though computationally expensive, captures these effects and is required for quantitative analysis.

**TIP3P:** A three-site water model (one oxygen, two hydrogens) with fixed charges and geometry. It is parametrised to work with AMBER force fields and gives a good balance of accuracy and speed. Its main weakness is slightly too-fast self-diffusion, but this is acceptable for most protein simulations.

### 2d. `genion` — neutralise and add salt

```bash
# First create a TPR (binary run input) from the current structure
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

# Replace water molecules with Na+ and Cl- ions
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
             -pname NA -nname CL -neutral -conc 0.15
```

**Why ions?** Proteins typically carry a net charge (e.g. 6D1Y is negatively charged). Simulating a charged system in a periodic box creates an infinite periodic array of net charge, which causes artefacts in electrostatic calculations. `genion` replaces random water molecules with Na⁺ or Cl⁻ ions to:

1. **Neutralise** the system (eliminate the net charge)
2. **Add 0.15 M NaCl** — approximately physiological salt concentration, which affects protein stability and surface electrostatics

The system is now complete: protein + water + ions, ready for simulation.

---

## Stage 3 — Energy Minimisation (`energy_minimize`)

```bash
gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

### What it does

The solvated system has many small structural problems: water molecules placed close to the protein may have slight overlaps, and the protein itself (reconstructed by PDBFixer) may have bond lengths or angles slightly off from ideal. These manifest as enormous repulsive forces between atoms.

Energy minimisation moves atoms to reduce the total potential energy of the system without propagating time. The **steepest descent** algorithm moves each atom a small step in the direction that reduces energy fastest, repeating until either:

- The maximum force on any atom drops below 1000 kJ/mol/nm (the `emtol` criterion), or
- 50,000 steps are exhausted

### Key MDP parameters

```
integrator  = steep      # steepest descent algorithm
emtol       = 1000.0     # convergence: max force < 1000 kJ/mol/nm
emstep      = 0.01       # initial step size in nm
nsteps      = 50000      # maximum steps
constraints = none       # no bond constraints during EM (uses all degrees of freedom)
```

No bond length constraints are applied during minimisation so the algorithm has full freedom to relieve strained geometries.

### Output

`em.gro` — minimised coordinates. `em.edr` contains the potential energy at each step; you can plot convergence with `gmx energy -f em.edr`. A successful minimisation shows the potential energy reaching a stable plateau and the max force falling below `emtol`.

---

## Stage 4 — Equilibration (`equilibrate`)

Equilibration brings the system to the correct temperature and pressure before the production run. It is split into two phases because coupling temperature and pressure simultaneously can be numerically unstable when starting from an energy-minimised structure with no initial velocities.

### 4a. NVT — constant volume, constant temperature

```
define         = -DPOSRES    # activate position restraints
integrator     = md          # velocity Verlet integrator
nsteps         = 50000       # 100 ps
tcoupl         = V-rescale   # velocity-rescaling thermostat
tc-grps        = Protein Non-Protein
ref_t          = 300 300     # target temperature for both groups
gen_vel        = yes         # generate initial velocities from Maxwell-Boltzmann
gen_temp       = 300
pcoupl         = no          # no pressure coupling yet
```

**Position restraints (`-DPOSRES`):** The heavy atoms of the protein are harmonically restrained to their minimised positions (force constant 1000 kJ/mol/nm²). This allows water and ions to move and relax freely around the protein while preventing the protein from drifting before the system has equilibrated thermally.

**V-rescale thermostat:** Velocities are periodically rescaled to maintain the target temperature. V-rescale adds a stochastic term that correctly samples the canonical (NVT) ensemble, unlike simple velocity rescaling which does not. Two coupling groups are used — Protein and Non-Protein (water + ions) — so each thermalises independently. This prevents the "hot solvent / cold solute" artefact where kinetic energy transferring slowly between groups leaves the protein artificially cold.

**Initial velocities:** Atoms begin with velocities drawn from the Maxwell-Boltzmann distribution at 300 K. The system immediately has thermal energy, but it needs time to equilibrate — the 100 ps NVT phase allows the kinetic and potential energy to reach a stable balance.

### 4b. NPT — constant pressure, constant temperature

```
define         = -DPOSRES    # protein still restrained
continuation   = yes         # continue from NVT checkpoint
gen_vel        = no          # reuse velocities from NVT
pcoupl         = Berendsen   # barostat: rescale box to reach target pressure
tau_p          = 2.0         # pressure coupling time constant (ps)
ref_p          = 1.0         # target pressure (bar)
compressibility = 4.5e-5     # water compressibility (bar⁻¹)
```

**Why add pressure coupling?** The NVT phase equilibrates temperature but the box volume is fixed. Real solution experiments happen at constant pressure, not constant volume. The NPT phase allows the box to expand or contract (by scaling coordinates) until the pressure stabilises at 1 bar.

**Berendsen barostat:** Chosen here for fast volume relaxation. It rescales the box and atomic coordinates at each step by a factor that drives pressure toward the target. Its weakness is that it does not rigorously sample the NPT ensemble (it suppresses volume fluctuations), but this is acceptable during equilibration where we only need the density to converge. The production run uses Parrinello-Rahman, which is thermodynamically correct.

**What to look for:** After NPT equilibration, the system density should be close to 1000 kg/m³ (the density of liquid water at 300 K and 1 bar). You can check with `gmx energy -f npt.edr` and selecting "Density".

---

## Stage 5 — Production Run (`production_run`)

```
integrator    = md
nsteps        = 500000        # 1 ns (500000 × 0.002 ps)
continuation  = yes
pcoupl        = Parrinello-Rahman   # correct NPT ensemble
constraints   = h-bonds             # LINCS algorithm for bonds to H
nstxout-compressed = 500            # save XTC frame every 500 steps = 1 ps
```

### The integrator and timestep

The **velocity Verlet** integrator (`md`) advances positions and velocities in discrete 2 fs steps:

```
v(t + dt/2) = v(t) + (F(t)/m) × (dt/2)
x(t + dt)   = x(t) + v(t + dt/2) × dt
v(t + dt)   = v(t + dt/2) + (F(t+dt)/m) × (dt/2)
```

The 2 fs timestep is chosen because the fastest motion in a protein is bond stretching (period ~10 fs). The Nyquist criterion requires at least 5 samples per period — in practice 2 fs is the largest step that maintains numerical stability for biomolecular systems. Going to 4 fs is possible but requires constraining more bond types and loses some accuracy.

### LINCS bond constraints

Bonds to hydrogen vibrate at ~3,300 cm⁻¹ (period ~10 fs). Integrating these with 2 fs steps is close to the stability limit. The **LINCS** algorithm (Linear Constraint Solver) freezes the length of all bonds involving hydrogen atoms, removing those degrees of freedom from the integration. This allows the 2 fs timestep while improving stability. The constraint does not significantly affect dynamics because hydrogen bond stretching is largely decoupled from biologically relevant motions.

### Parrinello-Rahman barostat

Unlike Berendsen, Parrinello-Rahman rigorously samples the isobaric-isothermal (NPT) ensemble, giving correct volume fluctuations. It is switched on only in production (not equilibration) because it can be numerically unstable if the system is far from equilibrium — which is why the Berendsen NPT equilibration phase is done first.

### Electrostatics — Particle Mesh Ewald (PME)

Long-range electrostatic interactions are computed with **PME**, which splits the Coulomb sum into a short-range real-space part (cut off at 1.0 nm) and a long-range reciprocal-space part computed on a grid via FFT. PME is O(N log N) in cost and gives essentially exact electrostatics — essential for charged proteins in ionic solution.

### Trajectory output

The compressed XTC format saves atom coordinates every `save_ps` picoseconds. A 10 ps save interval for a 1 ns run gives 100 frames — enough for basic analysis. For production-quality results:

- Use 10 ps saves for 10–100 ns runs
- Use 1–5 ps saves only if studying fast conformational changes (larger files)
- The raw XTC includes all atoms (protein + water + ions)
- `md_center.xtc` centres the protein and repairs PBC artefacts
- `md_fit.xtc` additionally removes overall rotation and translation — use this for all structural analysis

---

## Stage 6 — Trajectory Processing and Analysis

### Processing: centre and fit

```bash
# Step 1: Centre protein in box, make broken molecules whole
gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol -ur compact

# Step 2: Remove overall rotation and translation
gmx trjconv -s md.tpr -f md_center.xtc -o md_fit.xtc -fit rot+trans
```

**Why centre?** With PBC, the protein can drift and appear to "jump" across the box boundary between frames. `-center -pbc mol` ensures the protein stays in the middle of the box and all molecules are made whole.

**Why fit?** The protein undergoes diffusive rotation and translation (Brownian motion). These bulk motions carry no structural information and obscure internal fluctuations. Fitting each frame's backbone onto the reference structure (RMSD minimisation via the Kabsch algorithm) removes rotation and translation, so the remaining motion reflects internal protein dynamics only.

### RMSD — Root Mean Square Deviation

```bash
gmx rms -s md.tpr -f md_fit.xtc -o rmsd.xvg -tu ns
# Select Backbone for both fit group and RMSD group
```

**What it measures:** For each frame, the Cα backbone RMSD from the starting structure (after optimal superposition). Low RMSD = structure similar to start; high RMSD = structure has drifted.

**How to interpret:**

| RMSD pattern | Interpretation |
|---|---|
| Rises then plateaus at ~1–3 Å | Structure is stable; has equilibrated away from crystal contacts |
| Plateaus at < 1 Å | Very rigid protein or run too short to equilibrate |
| Continuously rising | Protein is unfolding, or simulation has errors — investigate |
| Multiple plateaus | Conformational transitions between distinct states |

An RMSD of 1.65 Å (as seen for 6D1Y at 1 ns) is normal — crystal structures are typically under lattice-induced strain and relax in solution.

### RMSF — Root Mean Square Fluctuation

```bash
gmx rmsf -s md.tpr -f md_fit.xtc -o rmsf.xvg -res
# Select C-alpha
```

**What it measures:** For each residue, the time-averaged Cα displacement from its mean position. Where RMSD tracks the whole protein over time, RMSF tracks each residue's flexibility over the whole trajectory.

High RMSF → flexible (loops, termini, disordered regions)
Low RMSF → rigid (buried core, secondary structure elements)

RMSF is directly comparable to crystallographic B-factors. The pipeline writes RMSF values into the B-factor column of `rmsf_bfactor.pdb` so PyMOL can colour the structure by dynamics with `spectrum b, blue_white_red`.

### Radius of Gyration (Rg)

```bash
gmx gyrate -s md.tpr -f md_fit.xtc -o gyrate.xvg
# Select Protein
```

**What it measures:** The mass-weighted RMS distance of all Cα atoms from the protein centre of mass — a measure of overall compactness.

A stable, globular protein shows a flat Rg trace. A rising Rg indicates expansion or unfolding. Rg also reveals breathing motions in multidomain proteins (domains moving apart and together).

### Energy and temperature

```bash
gmx energy -f md.edr -o energy.xvg
# Select Potential, Temperature, Pressure
```

Reads the binary energy file written during the simulation. The potential energy should fluctuate around a stable mean (drift indicates non-equilibration). Temperature should average 300 K with small fluctuations (±5 K is typical). Large temperature spikes indicate blown-up atoms — a sign of simulation instability.

### Hydrogen bond count

```bash
gmx hbond -s md.tpr -f md_fit.xtc -num hbond.xvg
# Select Protein, Protein
```

Counts intra-protein hydrogen bonds (donor–acceptor distance < 3.5 Å, D-H-A angle > 150°) at each frame. This is a sensitive indicator of secondary structure stability — a helix has characteristic NH→CO hydrogen bonds, and if the count drops sharply, a secondary structure element is breaking.

### DSSP — secondary structure per residue per frame

Run via MDAnalysis after the simulation:

```python
from MDAnalysis.analysis.dssp import DSSP
u = mda.Universe("md.tpr", "md_fit.xtc")
ds = DSSP(u.select_atoms("protein"))
ds.run()
# ds.results.dssp: shape (n_frames, n_residues), DSSP codes H/E/C/T/...
```

**DSSP codes:**

| Code | Structure |
|---|---|
| H | α-helix |
| G | 3₁₀-helix |
| I | π-helix |
| E | β-strand |
| B | isolated β-bridge |
| T | turn |
| S | bend |
| C | coil (none of the above) |

The output is a stacked area chart showing what fraction of residues are in each class at each time point. A stable protein maintains nearly constant secondary structure content throughout the run.

---

## Stage 7 — Visualisation

### Matplotlib plots

Six plots are generated and saved to `plots/`:

| Plot | What to look for |
|---|---|
| `rmsd.png` | Plateau = equilibrated; should stabilise within first 20–30% of run |
| `rmsf.png` | Peaks = flexible loops/termini; match to known functional regions |
| `rg.png` | Flat = compact & stable; drifting up = expanding or unfolding |
| `energy.png` | Flat mean with small fluctuations; temperature hovering at 300 K |
| `hbonds.png` | Roughly constant count; sharp drops = secondary structure disruption |
| `dssp.png` | Stable colored bands; shifts = secondary structure changes |

### PyMOL session

The `.pml` script loads `md_ref.pdb` as the base structure and appends all trajectory frames from `md_fit.xtc` as additional **states** on the same PyMOL object:

```python
load md_ref.pdb, protein          # state 1
load_traj md_fit.xtc, protein, 1  # states 2, 3, 4, ... N+1
```

`mset 1 -N` configures the movie to play states 1 through N in sequence (not `mset 1 xN` which would repeat state 1 N times with no motion).

Colouring by secondary structure (red helix, yellow strand, cyan loop) updates dynamically as you play through frames — you can watch helices form or break in real time.

The `rmsf_bfactor.pdb` file allows RMSF-based colouring: uncomment the lines in the `.pml` to load it, transfer B-factors to the trajectory object, and colour `blue_white_red` by B-factor. Blue = rigid, red = flexible.

---

## Running the Pipeline on a New PDB File

### Prerequisites

Install all dependencies once before your first run:

```bash
# 1. GROMACS (pick one)
brew install gromacs              # macOS
sudo apt install gromacs          # Ubuntu/Debian
conda install -c conda-forge gromacs
module load gromacs               # HPC — check your cluster's module name

# 2. Python packages (inside your venv)
source venv/bin/activate
pip install pdbfixer openmm       # required — fixes missing atoms
pip install MDAnalysis            # recommended — enables DSSP
pip install numpy matplotlib      # should already be installed
```

Verify everything is ready:

```bash
gmx --version            # should print GROMACS version
python3 -c "import pdbfixer; print('pdbfixer OK')"
python3 -c "import MDAnalysis; print('MDAnalysis OK')"
```

### Step 1 — Get your PDB file

```bash
# Download from RCSB (uses download_pdb.py)
python3 download_pdb.py --pdb 1ABC

# Or use an existing file
cp /path/to/my_structure.pdb Files/
```

**Before running MD, check the structure:**

```bash
python3 summarize_structures.py Files/1ABC.pdb
python3 summarize_structures.py Files/1ABC.pdb --missing
```

Look at:
- **Chains** — the pipeline simulates all chains together; if you want only one chain, extract it first
- **Missing residues** — PDBFixer will fill in missing *atoms* but cannot rebuild entire missing *loops*; large missing segments may need homology modelling first
- **Ligands** — the pipeline removes all HETATM records; if your biological question involves a ligand, you will need to parameterise it separately before using the pipeline
- **Resolution** — structures at < 2.5 Å are generally suitable; lower resolution structures may have more missing atoms

### Step 2 — Test run (confirm the pipeline works)

Always do a short test run before committing to a long simulation:

```bash
source venv/bin/activate
python3 run_md.py Files/1ABC.pdb \
    --ns 0.1 \
    --nvt-ps 10 \
    --npt-ps 10 \
    --ncores 4
```

This runs 10 ps NVT + 10 ps NPT + 100 ps production in a few minutes. Check that:
- All 7 stages complete without errors
- The `plots/` directory contains PNG files
- The `.pml` script exists and loads in PyMOL

If `pdb2gmx` fails with a missing atom error, confirm `pdbfixer` is installed and re-run (the prepare stage will be re-run automatically if `prepared.pdb` does not exist).

### Step 3 — Production run

Once the test run completes, run the full simulation. The pipeline is **resume-safe** — all completed stages are skipped automatically, so you can just re-run with the longer time:

```bash
python3 run_md.py Files/1ABC.pdb \
    --ns 50 \
    --ncores 4
```

**Recommended lengths by goal:**

| Goal | `--ns` |
|---|---|
| Pipeline test / stability check | 0.1–1 |
| Initial stability + flexibility | 10–50 |
| Side-chain rotamer sampling | 100–500 |
| Domain motion / conformational change | 100–1000 |

**To run on an HPC cluster**, write a job script that loads the GROMACS module, activates the venv, and runs the command. The `--ncores` flag controls OpenMP threads per MPI rank. For MPI-parallel GROMACS, consult your cluster documentation for the correct `mpirun` / `srun` invocation.

### Step 4 — Examine the outputs

After the run, all outputs are in `Files/1ABC_md/`:

```
Files/1ABC_md/
├── plots/
│   ├── 1ABC_summary.png      ← start here: all analyses in one figure
│   ├── 1ABC_rmsd.png
│   ├── 1ABC_rmsf.png
│   ├── 1ABC_rg.png
│   ├── 1ABC_energy.png
│   ├── 1ABC_hbonds.png
│   └── 1ABC_dssp.png         ← only if MDAnalysis installed
├── 1ABC_session.pml           ← open in PyMOL
├── rmsf_bfactor.pdb           ← colour structure by flexibility
├── md_fit.xtc                 ← fitted trajectory for further analysis
└── md.tpr                     ← binary topology (needed for gmx analysis tools)
```

**Check these first:**

1. Open `plots/1ABC_summary.png` — does the RMSD plateau? Is the temperature stable at 300 K?
2. Open PyMOL: `cd Files/1ABC_md && pymol 1ABC_session.pml`, then type `play`
3. Look at the RMSF plot: which residues are most flexible? Do they match known functional regions (active site loops, hinge regions)?

### Step 5 — Re-run analysis or visualisation only

If you want to change plot formatting or re-run just the analysis without re-running the simulation:

```bash
# Re-run analysis and visualisation only (simulation steps are all skipped)
python3 run_md.py Files/1ABC.pdb --steps analyze visualize
```

You can also run any GROMACS analysis tool manually on the saved trajectory and TPR:

```bash
cd Files/1ABC_md

# Solvent-accessible surface area over time
gmx sasa -s md.tpr -f md_fit.xtc -o sasa.xvg -surface "Protein" -output "Protein"

# Secondary structure per residue per frame (if gmx dssp is available, GROMACS ≥ 2023)
gmx dssp -s md.tpr -f md_fit.xtc -o dssp.dat

# Distance between two residue groups over time
gmx distance -s md.tpr -f md_fit.xtc -oall dist.xvg -select "com of resid 50 plus com of resid 120"

# Cluster trajectory into representative conformations
gmx cluster -s md.tpr -f md_fit.xtc -method gromos -cutoff 0.2 -cl clusters.pdb
```

### Common issues and fixes

| Error | Cause | Fix |
|---|---|---|
| `atom CG not found in residue LYS` | Missing side-chain atoms in input PDB | `pip install pdbfixer openmm` (if not done) |
| `pdb2gmx` fails with unknown residue | Non-standard residue (modified AA, ligand) | Remove the residue from the PDB or parameterise it separately |
| `mdrun` crashes immediately | System has bad contacts after solvation | The minimisation should handle this; if it persists, try increasing `emtol` to 10000 or reducing `emstep` |
| RMSD drifts continuously upward | Protein unfolding or simulation instability | Check that equilibration converged; try longer equilibration (`--nvt-ps 200 --npt-ps 200`) |
| Temperature fluctuations > ±20 K | Timestep too large or constraints failing | Reduce `--save-ps` is not the fix; check the MDP logs in the `_md/` directory |
| PyMOL plays with no motion | Wrong `mset` syntax or trajectory not loaded | Ensure the `.pml` uses `mset 1 -N` not `mset 1 xN`; confirm `md_fit.xtc` exists |
| Very slow simulation | Too many atoms or too few cores | Check system size with `gmx editconf -f solv_ions.gro -o /dev/null`; increase `--ncores` |
