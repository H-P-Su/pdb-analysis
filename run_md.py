#!/usr/bin/env python3
"""
run_md.py — GROMACS molecular dynamics pipeline for protein structures.

Pipeline stages (run all by default, or select with --steps):
  prepare   Clean PDB: strip HETATM/ANISOU; add missing atoms via PDBFixer if available
  setup     Build topology (pdb2gmx), define box (editconf), solvate, add ions (genion)
  minimize  Steepest-descent energy minimisation (~50,000 steps)
  equil     NVT equilibration (heat to target T) → NPT equilibration (equilibrate density)
  run       Production MD (default 10 ns; save XTC every 10 ps)
  analyze   RMSD, RMSF, radius of gyration, energy/temperature, H-bond count, DSSP
  visualize Matplotlib time-series plots + PyMOL session script (.pml)

Requirements — install into venv before use:
  GROMACS ≥ 2019  (LGPL v2.1)  https://www.gromacs.org
                                macOS:  brew install gromacs
                                Ubuntu: sudo apt install gromacs
                                HPC:    module load gromacs
  numpy           (BSD)         pip install numpy
  matplotlib      (BSD)         pip install matplotlib
  MDAnalysis      (GPL v2)      pip install MDAnalysis        [recommended — enables DSSP]
  pdbfixer        (MIT)         pip install pdbfixer openmm   [optional  — better structure prep]

Usage:
  python3 run_md.py protein.pdb                            # full pipeline, 10 ns
  python3 run_md.py protein.pdb --ns 100                   # 100 ns production run
  python3 run_md.py protein.pdb --steps prepare setup      # only build inputs
  python3 run_md.py protein.pdb --steps analyze visualize  # post-process existing run
  python3 run_md.py protein.pdb --ff charmm36m-iua --water tip4p
  python3 run_md.py protein.pdb --ncores 8 --gpu           # parallel / GPU run

All outputs are written to <pdb_stem>_md/ (one subdirectory per structure).
Re-running the script skips completed steps automatically (resume-safe).
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend — works on headless servers
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

log = logging.getLogger("run_md")

# ─── constants ────────────────────────────────────────────────────────────────

SUPPORTED_FORCEFIELDS = {
    "amber99sb-ildn": "AMBER99SB-ILDN (recommended for globular proteins)",
    "amber03":        "AMBER03",
    "charmm36m-iua":  "CHARMM36m (recommended for IDPs and membrane proteins)",
    "oplsaa":         "OPLS-AA/L",
    "gromos54a7":     "GROMOS54a7",
}
WATER_MODELS = {
    "tip3p":  "TIP3P  — default; works with AMBER/CHARMM force fields",
    "tip4p":  "TIP4P  — 4-site; better bulk water properties",
    "spc":    "SPC",
    "spce":   "SPC/E  — better diffusion; use with GROMOS",
}

DEFAULT_FF     = "amber99sb-ildn"
DEFAULT_WATER  = "tip3p"
DEFAULT_NS     = 10       # production run length in ns
DEFAULT_DT     = 0.002    # integration timestep in ps
DEFAULT_BOX    = 1.0      # min protein-to-box-edge distance in nm
DEFAULT_CONC   = 0.15     # NaCl concentration in mol/L
DEFAULT_TEMP   = 300.0    # simulation temperature in K
DEFAULT_SAVEPS = 10.0     # trajectory save interval in ps

ALL_STEPS = ["prepare", "setup", "minimize", "equil", "run", "analyze", "visualize"]

# ─── OpenMM constants ─────────────────────────────────────────────────────────

OPENMM_FORCEFIELDS = {
    "amber14-all":   "AMBER14 all-atom / ff14SB (recommended for OpenMM)",
    "amber99sbildn": "AMBER99SB-ILDN — equivalent to GROMACS amber99sb-ildn",
    "charmm36":      "CHARMM36m — IDPs, nucleic acids, membrane proteins",
}
OPENMM_WATER_MODELS = {
    "tip3p":   "TIP3P  — standard; compatible with amber14 and amber99sbildn",
    "tip4pew": "TIP4P-Ew — better bulk water; use with amber14",
    "spce":    "SPC/E  — use with charmm36",
}
# Maps (ff, water) → list of OpenMM ForceField XML file names
_OMM_FF_XMLS: dict = {
    ("amber14-all",   "tip3p"):   ["amber14-all.xml", "amber14/tip3p.xml"],
    ("amber14-all",   "tip4pew"): ["amber14-all.xml", "amber14/tip4pew.xml"],
    ("amber14-all",   "spce"):    ["amber14-all.xml", "amber14/spce.xml"],
    ("amber99sbildn", "tip3p"):   ["amber99sbildn.xml", "tip3p.xml"],
    ("amber99sbildn", "tip4pew"): ["amber99sbildn.xml", "tip4pew.xml"],
    ("charmm36",      "tip3p"):   ["charmm36.xml", "charmm36/water.xml"],
    ("charmm36",      "spce"):    ["charmm36.xml", "charmm36/water.xml"],
}
DEFAULT_ENGINE    = "gromacs"
DEFAULT_OMM_FF    = "amber14-all"
DEFAULT_OMM_WATER = "tip3p"

# ─── MDP file templates ───────────────────────────────────────────────────────
# These are written to disk before each stage. Format placeholders are filled
# by _write_mdp() from MDConfig properties.

MDP_MINIM = """\
; Energy minimisation — steepest descent until F_max < 1000 kJ/mol/nm
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000

; Output
nstenergy       = 500
nstlog          = 500
nstxout         = 0
nstvout         = 0

; Neighbour search
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 1
rcoulomb        = 1.0
rvdw            = 1.0

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16

; Van der Waals
DispCorr        = EnerPres

; No constraints during minimisation
constraints     = none
"""

MDP_NVT = """\
; NVT equilibration — {ns_nvt} ps at {temp} K (constant volume)
; Position restraints on heavy protein atoms keep the protein fixed
; while water and ions equilibrate around it.
define          = -DPOSRES
integrator      = md
dt              = {dt}
nsteps          = {steps_nvt}

; Output
nstxout         = 500
nstvout         = 500
nstenergy       = 500
nstlog          = 500

; Bonds
continuation    = no
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbour search
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16

; Van der Waals
DispCorr        = EnerPres

; Temperature coupling (V-rescale = velocity rescaling with stochastic term)
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1 0.1
ref_t           = {temp} {temp}

; No pressure coupling in NVT
pcoupl          = no

; PBC
pbc             = xyz

; Generate initial velocities from Maxwell-Boltzmann distribution
gen_vel         = yes
gen_temp        = {temp}
gen_seed        = -1
"""

MDP_NPT = """\
; NPT equilibration — {ns_npt} ps at {temp} K, 1 bar (constant pressure)
; Position restraints remain on protein heavy atoms.
; Berendsen barostat used here for fast volume relaxation.
define          = -DPOSRES
integrator      = md
dt              = {dt}
nsteps          = {steps_npt}

; Output
nstxout         = 500
nstvout         = 500
nstenergy       = 500
nstlog          = 500

; Bonds
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbour search
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16

; Van der Waals
DispCorr        = EnerPres

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1 0.1
ref_t           = {temp} {temp}

; Pressure coupling (Berendsen for fast equilibration)
pcoupl          = Berendsen
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5
refcoord_scaling = com

; PBC
pbc             = xyz

; Continue from NVT — do not regenerate velocities
gen_vel         = no
"""

MDP_MD = """\
; Production MD — {ns_prod} ns at {temp} K, 1 bar
; Parrinello-Rahman barostat samples true NPT ensemble.
; No position restraints.
integrator      = md
dt              = {dt}
nsteps          = {steps_prod}

; Save trajectory every {save_ps} ps
nstxout             = {save_steps}
nstvout             = {save_steps}
nstenergy           = {save_steps}
nstlog              = {save_steps}
nstxout-compressed  = {save_steps}
compressed-x-grps   = System

; Bonds
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbour search
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16

; Van der Waals
DispCorr        = EnerPres

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1 0.1
ref_t           = {temp} {temp}

; Pressure coupling (Parrinello-Rahman for correct NPT ensemble)
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5

; PBC
pbc             = xyz

; Continue from NPT
gen_vel         = no
"""

# ─── configuration ────────────────────────────────────────────────────────────

@dataclass
class MDConfig:
    """All parameters and file paths for one MD pipeline run."""
    pdb_in:   Path
    workdir:  Path
    ff:       str   = DEFAULT_FF
    water:    str   = DEFAULT_WATER
    ns_prod:  float = DEFAULT_NS
    dt:       float = DEFAULT_DT
    temp:     float = DEFAULT_TEMP
    box_dist: float = DEFAULT_BOX
    ion_conc: float = DEFAULT_CONC
    ns_nvt:   float = 0.1    # NVT equilibration in ns
    ns_npt:   float = 0.1    # NPT equilibration in ns
    save_ps:  float = DEFAULT_SAVEPS
    ncores:   int   = 0      # 0 = let GROMACS/OpenMM auto-detect
    gpu:      bool  = False
    steps:    list  = field(default_factory=lambda: list(ALL_STEPS))
    engine:   str   = DEFAULT_ENGINE

    # All output paths derived from workdir
    prepared:  Path = field(init=False)
    processed: Path = field(init=False)
    topol:     Path = field(init=False)
    newbox:    Path = field(init=False)
    solv:      Path = field(init=False)
    solv_ions: Path = field(init=False)
    em_tpr:    Path = field(init=False)
    em_gro:    Path = field(init=False)
    nvt_tpr:   Path = field(init=False)
    nvt_gro:   Path = field(init=False)
    npt_tpr:   Path = field(init=False)
    npt_gro:   Path = field(init=False)
    md_tpr:    Path = field(init=False)
    md_gro:    Path = field(init=False)
    md_xtc:    Path = field(init=False)
    md_edr:    Path = field(init=False)
    md_center: Path = field(init=False)
    md_fit:    Path = field(init=False)
    plots_dir: Path = field(init=False)

    def __post_init__(self):
        d = self.workdir
        self.prepared  = d / "prepared.pdb"
        self.processed = d / "processed.gro"
        self.topol     = d / "topol.top"
        self.newbox    = d / "newbox.gro"
        self.solv      = d / "solv.gro"
        self.solv_ions = d / "solv_ions.gro"
        self.em_tpr    = d / "em.tpr"
        self.em_gro    = d / "em.gro"
        self.nvt_tpr   = d / "nvt.tpr"
        self.nvt_gro   = d / "nvt.gro"
        self.npt_tpr   = d / "npt.tpr"
        self.npt_gro   = d / "npt.gro"
        self.md_tpr    = d / "md.tpr"
        self.md_gro    = d / "md.gro"
        self.md_xtc    = d / "md.xtc"
        self.md_edr    = d / "md.edr"
        self.md_center = d / "md_center.xtc"
        self.md_fit    = d / "md_fit.xtc"
        self.plots_dir = d / "plots"

        # OpenMM-specific output paths
        self.omm_topology:  Path = d / "omm_topology.pdb"   # solvated system topology
        self.omm_system:    Path = d / "omm_system.xml"      # serialised OpenMM System
        self.omm_em_pdb:    Path = d / "em.pdb"              # post-minimisation coords
        self.omm_nvt_pdb:   Path = d / "nvt.pdb"             # post-NVT coords
        self.omm_npt_pdb:   Path = d / "npt.pdb"             # post-NPT coords
        self.omm_md_dcd:    Path = d / "md_traj.dcd"         # raw production trajectory
        self.omm_md_fit:    Path = d / "md_fit.dcd"          # backbone-fitted protein-only DCD
        self.omm_energy_csv:Path = d / "omm_energy.csv"      # StateDataReporter output

    # Derived step counts from ns and dt
    @property
    def steps_nvt(self) -> int:
        return round(self.ns_nvt * 1000 / self.dt)

    @property
    def steps_npt(self) -> int:
        return round(self.ns_npt * 1000 / self.dt)

    @property
    def steps_prod(self) -> int:
        return round(self.ns_prod * 1000 / self.dt)

    @property
    def save_steps(self) -> int:
        return max(1, round(self.save_ps / self.dt))

    @property
    def n_frames(self) -> int:
        return max(1, round(self.ns_prod * 1000 / self.save_ps))

    @property
    def active_trajectory(self) -> Path:
        """Fitted trajectory file (engine-agnostic), for PyMOL and analysis."""
        if self.engine == "openmm":
            return self.omm_md_fit if self.omm_md_fit.exists() else self.omm_md_dcd
        return self.md_fit if self.md_fit.exists() else self.md_xtc

    @property
    def mdrun_flags(self) -> list:
        flags = []
        if self.ncores > 0:
            flags += ["-ntmpi", "1", "-ntomp", str(self.ncores)]
        if self.gpu:
            flags += ["-gpu_id", "0"]
        return flags


# ─── GROMACS runner ───────────────────────────────────────────────────────────

GMX: str = ""   # filled by _find_gmx() in main()


def _find_gmx() -> str:
    """Return the GROMACS executable name, or exit with an install hint."""
    for name in ("gmx", "gmx_mpi", "gmx_d", "gmx_mpi_d"):
        if shutil.which(name):
            return name
    sys.exit(
        "\n[error] GROMACS not found in PATH.\n\n"
        "Install GROMACS and add it to PATH before running this script:\n"
        "  macOS (Homebrew):  brew install gromacs\n"
        "  Ubuntu/Debian:     sudo apt install gromacs\n"
        "  Conda:             conda install -c conda-forge gromacs\n"
        "  HPC cluster:       module load gromacs\n"
        "  From source:       https://manual.gromacs.org/current/install-guide\n"
    )


def gmx(*args, stdin: str = "", check: bool = True,
         cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    """Run a GROMACS subcommand, log it, and raise RuntimeError on failure."""
    cmd = [GMX] + [str(a) for a in args]
    log.info("  $ %s", " ".join(str(x) for x in cmd))
    if stdin:
        log.debug("    stdin: %r", stdin)
    result = subprocess.run(
        cmd,
        input=stdin.encode(),
        capture_output=True,
        cwd=str(cwd) if cwd else None,
    )
    if check and result.returncode != 0:
        stderr = result.stderr.decode(errors="replace")
        stdout = result.stdout.decode(errors="replace")
        log.error(
            "GROMACS command failed (exit %d):\n  %s\n--- stderr ---\n%s",
            result.returncode, " ".join(cmd), stderr[-4000:]
        )
        if stdout.strip():
            log.debug("--- stdout ---\n%s", stdout[-2000:])
        raise RuntimeError(f"gmx {args[0]} failed (exit {result.returncode})")
    return result


# ─── utility helpers ──────────────────────────────────────────────────────────

def _write_mdp(path: Path, template: str, **kwargs) -> None:
    path.write_text(template.format(**kwargs))
    log.debug("  Wrote %s", path.name)


def _skip(path: Path, label: str) -> bool:
    """Return True (and log) if path already exists with content."""
    if path.exists() and path.stat().st_size > 0:
        log.info("  [skip] %s already exists (%s)", path.name, label)
        return True
    return False


def parse_xvg(path: Path) -> tuple[np.ndarray, list[str]]:
    """
    Parse a GROMACS XVG file into a NumPy array.

    Returns (data, column_labels) where:
      data   — shape (n_frames, n_columns); column 0 is always time
      labels — list of y-axis legend strings from '@ sN legend' lines
    """
    rows, labels = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("@"):
                if "legend" in line:
                    parts = line.split('"')
                    if len(parts) >= 2:
                        labels.append(parts[1])
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    return np.array(rows) if rows else np.empty((0, 2)), labels


# ─── stage 1: prepare ─────────────────────────────────────────────────────────

def prepare_structure(cfg: MDConfig) -> None:
    """
    Stage 1 — Clean and fix the input PDB.

    With PDBFixer (pip install pdbfixer openmm):
      - Fills missing heavy atoms and short loops
      - Replaces non-standard residues with standard equivalents
      - Removes all HETATM records (water, ligands, ions)

    Without PDBFixer (fallback):
      - Strips HETATM, ANISOU, CONECT, and REMARK records
      - Keeps only ATOM, TER, MODEL, ENDMDL, CRYST1 lines

    Note: ligands are removed in both cases because standard force fields
    (AMBER, CHARMM, OPLS) do not have parameters for arbitrary small molecules.
    To include a ligand, generate its topology separately with ACPYPE or ATB
    and merge it manually before running --steps setup.
    """
    log.info("[prepare] Cleaning PDB → %s", cfg.prepared.name)
    if _skip(cfg.prepared, "prepare"):
        return

    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
        log.info("  PDBFixer detected — fixing missing atoms and residues")
        fixer = PDBFixer(str(cfg.pdb_in))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        with open(cfg.prepared, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh)
        log.info("  PDBFixer: done")
    except ImportError:
        log.info("  pdbfixer not installed — using basic HETATM strip")
        _basic_clean_pdb(cfg.pdb_in, cfg.prepared)


def _basic_clean_pdb(src: Path, dst: Path) -> None:
    keep = ("ATOM", "TER", "MODEL", "ENDMDL", "CRYST1")
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if any(line.startswith(p) for p in keep):
                fo.write(line)
        fo.write("END\n")
    log.info("  Stripped HETATM / ANISOU records")


# ─── stage 2: setup ───────────────────────────────────────────────────────────

def setup_topology(cfg: MDConfig) -> None:
    """
    Stage 2 — Build the simulation system from the cleaned PDB.

    Sub-steps (each skipped if output already exists):
      pdb2gmx  — assigns force-field atom types, adds hydrogen atoms,
                 writes topol.top and posre.itp (position restraint file)
      editconf — places protein in a dodecahedral periodic box with
                 ≥ cfg.box_dist nm between protein and box edge
      solvate  — fills the box with water (SPC/E or TIP3P geometry)
      genion   — replaces water molecules with Na+ and Cl- to neutralise
                 the system and reach the target salt concentration

    GROMACS atom group numbers used (standard for protein-in-water systems):
      0 = System        1 = Protein      3 = C-alpha
      4 = Backbone      SOL = solvent     (passed by name to genion)
    """
    log.info("[setup] Building solvated topology")

    # pdb2gmx: force field → topology + initial coordinates
    if not _skip(cfg.processed, "pdb2gmx"):
        log.info("  pdb2gmx (%s / %s water) …", cfg.ff, cfg.water)
        gmx(
            "pdb2gmx",
            "-f", cfg.prepared,
            "-o", cfg.processed,
            "-p", cfg.topol,
            "-i", cfg.workdir / "posre.itp",
            "-ff", cfg.ff,
            "-water", cfg.water,
            "-ignh",          # regenerate all H atoms from scratch
            cwd=cfg.workdir,
        )

    # editconf: define periodic box
    if not _skip(cfg.newbox, "editconf"):
        log.info("  editconf (dodecahedron, d=%.1f nm) …", cfg.box_dist)
        gmx(
            "editconf",
            "-f", cfg.processed,
            "-o", cfg.newbox,
            "-c",               # centre protein
            "-d", str(cfg.box_dist),
            "-bt", "dodecahedron",
            cwd=cfg.workdir,
        )

    # solvate: add water
    if not _skip(cfg.solv, "solvate"):
        log.info("  solvate …")
        gmx(
            "solvate",
            "-cp", cfg.newbox,
            "-cs", "spc216.gro",
            "-o", cfg.solv,
            "-p", cfg.topol,
            cwd=cfg.workdir,
        )

    # genion: add Na+/Cl- ions
    if not _skip(cfg.solv_ions, "genion"):
        log.info("  genion (%.2f M NaCl, neutralise) …", cfg.ion_conc)
        ions_mdp = cfg.workdir / "ions.mdp"
        ions_tpr = cfg.workdir / "ions.tpr"
        _write_mdp(ions_mdp, MDP_MINIM)
        gmx("grompp",
            "-f", ions_mdp, "-c", cfg.solv, "-p", cfg.topol,
            "-o", ions_tpr, "-maxwarn", "2",
            cwd=cfg.workdir)
        gmx("genion",
            "-s", ions_tpr, "-o", cfg.solv_ions, "-p", cfg.topol,
            "-pname", "NA", "-nname", "CL", "-neutral", "-conc", str(cfg.ion_conc),
            stdin="SOL\n",     # replace water molecules with ions
            cwd=cfg.workdir)


# ─── stage 3: minimize ────────────────────────────────────────────────────────

def energy_minimize(cfg: MDConfig) -> None:
    """
    Stage 3 — Steepest-descent energy minimisation.

    Removes steric clashes introduced during solvation.
    Runs up to 50,000 steps or until max force < 1000 kJ/mol/nm.
    """
    log.info("[minimize] Energy minimisation")
    minim_mdp = cfg.workdir / "minim.mdp"
    _write_mdp(minim_mdp, MDP_MINIM)

    if not _skip(cfg.em_tpr, "grompp-em"):
        gmx("grompp",
            "-f", minim_mdp, "-c", cfg.solv_ions, "-p", cfg.topol,
            "-o", cfg.em_tpr, cwd=cfg.workdir)

    if not _skip(cfg.em_gro, "mdrun-em"):
        log.info("  mdrun (minimisation) …")
        gmx("mdrun", "-v", "-deffnm", "em", *cfg.mdrun_flags, cwd=cfg.workdir)


# ─── stage 4: equil ───────────────────────────────────────────────────────────

def equilibrate(cfg: MDConfig) -> None:
    """
    Stage 4 — Two-phase equilibration with heavy-atom position restraints.

    NVT (constant N, V, T):
      Protein is restrained; water and ions move freely. System is heated
      from 0 K to cfg.temp over cfg.ns_nvt ps using V-rescale thermostat.

    NPT (constant N, P, T):
      Pressure coupling (Berendsen) is turned on. Box volume relaxes to
      equilibrium density (~1000 kg/m³ for TIP3P water at 300 K).
    """
    log.info("[equil] NVT + NPT equilibration")

    # ── NVT ───────────────────────────────────────────────────────────────────
    nvt_mdp = cfg.workdir / "nvt.mdp"
    _write_mdp(nvt_mdp, MDP_NVT,
               ns_nvt=int(cfg.ns_nvt * 1000), dt=cfg.dt,
               steps_nvt=cfg.steps_nvt, temp=int(cfg.temp))

    if not _skip(cfg.nvt_tpr, "grompp-nvt"):
        gmx("grompp",
            "-f", nvt_mdp, "-c", cfg.em_gro, "-r", cfg.em_gro,
            "-p", cfg.topol, "-o", cfg.nvt_tpr, "-maxwarn", "1",
            cwd=cfg.workdir)

    if not _skip(cfg.nvt_gro, "mdrun-nvt"):
        log.info("  NVT mdrun (%d ps) …", int(cfg.ns_nvt * 1000))
        gmx("mdrun", "-v", "-deffnm", "nvt", *cfg.mdrun_flags, cwd=cfg.workdir)

    # ── NPT ───────────────────────────────────────────────────────────────────
    npt_mdp = cfg.workdir / "npt.mdp"
    _write_mdp(npt_mdp, MDP_NPT,
               ns_npt=int(cfg.ns_npt * 1000), dt=cfg.dt,
               steps_npt=cfg.steps_npt, temp=int(cfg.temp))

    if not _skip(cfg.npt_tpr, "grompp-npt"):
        gmx("grompp",
            "-f", npt_mdp, "-c", cfg.nvt_gro, "-r", cfg.nvt_gro,
            "-t", cfg.workdir / "nvt.cpt",
            "-p", cfg.topol, "-o", cfg.npt_tpr, "-maxwarn", "1",
            cwd=cfg.workdir)

    if not _skip(cfg.npt_gro, "mdrun-npt"):
        log.info("  NPT mdrun (%d ps) …", int(cfg.ns_npt * 1000))
        gmx("mdrun", "-v", "-deffnm", "npt", *cfg.mdrun_flags, cwd=cfg.workdir)


# ─── stage 5: production run ──────────────────────────────────────────────────

def production_run(cfg: MDConfig) -> None:
    """
    Stage 5 — Production MD.

    No position restraints. Full Parrinello-Rahman barostat (correct NPT
    ensemble). Compressed XTC trajectory saved every cfg.save_ps picoseconds.
    """
    log.info("[run] Production MD (%.1f ns, %d steps)", cfg.ns_prod, cfg.steps_prod)
    md_mdp = cfg.workdir / "md.mdp"
    _write_mdp(md_mdp, MDP_MD,
               ns_prod=cfg.ns_prod, dt=cfg.dt,
               steps_prod=cfg.steps_prod, temp=int(cfg.temp),
               save_ps=cfg.save_ps, save_steps=cfg.save_steps)

    if not _skip(cfg.md_tpr, "grompp-md"):
        gmx("grompp",
            "-f", md_mdp, "-c", cfg.npt_gro, "-t", cfg.workdir / "npt.cpt",
            "-p", cfg.topol, "-o", cfg.md_tpr,
            cwd=cfg.workdir)

    if not _skip(cfg.md_xtc, "mdrun-md"):
        log.info("  mdrun (%.1f ns, dt=%.3f ps) …", cfg.ns_prod, cfg.dt)
        gmx("mdrun", "-v", "-deffnm", "md", *cfg.mdrun_flags, cwd=cfg.workdir)


# ─── OpenMM helpers ───────────────────────────────────────────────────────────

def _omm_platform(cfg: MDConfig):
    """Return the fastest available OpenMM Platform given --gpu / --ncores."""
    import openmm
    if cfg.gpu:
        for name in ("CUDA", "OpenCL"):
            try:
                return openmm.Platform.getPlatformByName(name)
            except Exception:
                pass
        log.warning("  No GPU platform found — falling back to CPU")
    plat = openmm.Platform.getPlatformByName("CPU")
    if cfg.ncores > 0:
        plat.setPropertyDefaultValue("Threads", str(cfg.ncores))
    return plat


def _omm_check_import() -> None:
    """Exit with a helpful message if OpenMM is not installed."""
    try:
        import openmm  # noqa: F401
    except ImportError:
        sys.exit(
            "\n[error] OpenMM not installed.\n\n"
            "Install with:\n"
            "  pip install openmm pdbfixer\n"
            "Verify:\n"
            "  python -c 'import openmm; print(openmm.__version__)'\n"
        )


def _omm_load_pdb(path: Path):
    """Load an OpenMM PDBFile and return (topology, positions)."""
    from openmm.app import PDBFile
    pdb = PDBFile(str(path))
    return pdb.topology, pdb.positions


# ─── OpenMM stage 2: setup ────────────────────────────────────────────────────

def omm_setup(cfg: MDConfig) -> None:
    """
    OpenMM Stage 2 — Build solvated system.

    Steps:
      - Load the prepared PDB (from stage 1)
      - Add missing hydrogens at pH 7.0
      - Solvate in a cubic box (padding = cfg.box_dist nm each side)
      - Add Na+/Cl- ions to neutralise and reach cfg.ion_conc mol/L
      - Save solvated topology as omm_topology.pdb
      - Serialise the OpenMM System to omm_system.xml for later stages
    """
    _omm_check_import()
    log.info("[setup/openmm] Building solvated OpenMM system")
    if _skip(cfg.omm_topology, "omm-setup"):
        return

    import openmm
    import openmm.unit as unit
    from openmm.app import PDBFile, ForceField, Modeller, PME, HBonds

    key = (cfg.ff, cfg.water)
    if key not in _OMM_FF_XMLS:
        supported = ", ".join(f"{f}/{w}" for f, w in _OMM_FF_XMLS)
        sys.exit(
            f"\n[error] Unsupported force field / water model combination for OpenMM: "
            f"{cfg.ff} / {cfg.water}\n"
            f"Supported pairs: {supported}\n"
        )
    ff_xmls = _OMM_FF_XMLS[key]
    log.info("  Force field XMLs: %s", ", ".join(ff_xmls))

    forcefield = ForceField(*ff_xmls)
    pdb        = PDBFile(str(cfg.prepared))
    modeller   = Modeller(pdb.topology, pdb.positions)

    log.info("  Adding missing hydrogens (pH 7.0) …")
    modeller.addHydrogens(forcefield, pH=7.0)

    log.info("  Solvating (padding=%.1f nm, [NaCl]=%.2f M, neutralise) …",
             cfg.box_dist, cfg.ion_conc)
    modeller.addSolvent(
        forcefield,
        model=cfg.water,
        padding=cfg.box_dist * unit.nanometers,
        ionicStrength=cfg.ion_conc * unit.molar,
        neutralize=True,
    )
    n_atoms = modeller.topology.getNumAtoms()
    log.info("  Solvated system: %d atoms", n_atoms)

    # Save solvated topology (reference for MDAnalysis and PyMOL)
    with open(cfg.omm_topology, "w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh)
    log.info("  Topology → %s", cfg.omm_topology.name)

    # Build System and serialise it (avoid re-running ForceField.createSystem later)
    log.info("  Creating system (PME, HBonds constraints, 1 nm cutoff) …")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=HBonds,
    )
    with open(cfg.omm_system, "w") as fh:
        fh.write(openmm.XmlSerializer.serialize(system))
    log.info("  System → %s", cfg.omm_system.name)


# ─── OpenMM stage 3: minimize ─────────────────────────────────────────────────

def omm_minimize(cfg: MDConfig) -> None:
    """
    OpenMM Stage 3 — Energy minimisation with L-BFGS.

    Runs until convergence (max 50 000 iterations) or until the energy
    gradient is below the OpenMM default tolerance (10 kJ/mol/nm).
    Saves minimised coordinates to em.pdb.
    """
    _omm_check_import()
    log.info("[minimize/openmm] Energy minimisation (L-BFGS)")
    if _skip(cfg.omm_em_pdb, "omm-minimize"):
        return

    import openmm
    import openmm.unit as unit
    from openmm.app import PDBFile, Simulation, LocalEnergyMinimizer

    topology, positions = _omm_load_pdb(cfg.omm_topology)
    system = openmm.XmlSerializer.deserialize(cfg.omm_system.read_text())

    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temp * unit.kelvin, 1.0 / unit.picosecond,
        cfg.dt * unit.picoseconds)
    sim = Simulation(topology, system, integrator, _omm_platform(cfg))
    sim.context.setPositions(positions)

    e0 = (sim.context.getState(getEnergy=True)
          .getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
    log.info("  Initial energy: %.1f kJ/mol", e0)
    LocalEnergyMinimizer.minimize(sim.context, maxIterations=50000)
    e1 = (sim.context.getState(getEnergy=True)
          .getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
    log.info("  Final energy:   %.1f kJ/mol", e1)

    state = sim.context.getState(getPositions=True)
    with open(cfg.omm_em_pdb, "w") as fh:
        PDBFile.writeFile(topology, state.getPositions(), fh)
    log.info("  Minimised structure → %s", cfg.omm_em_pdb.name)


# ─── OpenMM stage 4: equilibrate ──────────────────────────────────────────────

def omm_equilibrate(cfg: MDConfig) -> None:
    """
    OpenMM Stage 4 — Two-phase equilibration with backbone position restraints.

    NVT (100 ps default):
      LangevinMiddleIntegrator at cfg.temp K.  Backbone heavy atoms are
      restrained with a harmonic force (k = 1000 kJ/mol/nm²) so water
      equilibrates around the fixed protein.

    NPT (100 ps default):
      MonteCarloBarostat added.  Volume relaxes to equilibrium density.
      Restraints remain active.  State saved to nvt.chk / npt.chk for
      resumption in stage 5.
    """
    _omm_check_import()
    log.info("[equil/openmm] NVT + NPT equilibration")
    if _skip(cfg.omm_npt_pdb, "omm-equil"):
        return

    import openmm
    import openmm.unit as unit
    from openmm.app import PDBFile, Simulation, StateDataReporter

    topology, _ = _omm_load_pdb(cfg.omm_topology)
    _, positions = _omm_load_pdb(cfg.omm_em_pdb)
    system = openmm.XmlSerializer.deserialize(cfg.omm_system.read_text())

    # Position restraints on protein backbone heavy atoms
    restraint = openmm.CustomExternalForce(
        "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter(
        "k", 1000.0 * unit.kilojoules_per_mole / unit.nanometers**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    n_restrained = 0
    for atom in topology.atoms():
        if (atom.residue.name not in ("HOH", "WAT", "NA", "CL", "Na+", "Cl-")
                and atom.element is not None
                and atom.element.symbol not in ("H", "D")):
            p = positions[atom.index]
            restraint.addParticle(atom.index, [p.x, p.y, p.z])
            n_restrained += 1
    system.addForce(restraint)
    log.info("  Position restraints on %d protein heavy atoms", n_restrained)

    # ── NVT ───────────────────────────────────────────────────────────────────
    nvt_chk = cfg.workdir / "nvt.chk"
    if not _skip(cfg.omm_nvt_pdb, "omm-nvt"):
        nvt_steps = cfg.steps_nvt
        log.info("  NVT (%d steps = %d ps) …", nvt_steps, int(cfg.ns_nvt * 1000))
        integrator = openmm.LangevinMiddleIntegrator(
            cfg.temp * unit.kelvin, 1.0 / unit.picosecond,
            cfg.dt * unit.picoseconds)
        integrator.setRandomNumberSeed(42)
        sim = Simulation(topology, system, integrator, _omm_platform(cfg))
        sim.context.setPositions(positions)
        sim.context.setVelocitiesToTemperature(cfg.temp * unit.kelvin, 42)
        sim.reporters.append(StateDataReporter(
            str(cfg.workdir / "nvt_log.csv"), max(500, nvt_steps // 100),
            step=True, time=True, temperature=True, potentialEnergy=True))
        sim.step(nvt_steps)
        state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
        with open(cfg.omm_nvt_pdb, "w") as fh:
            PDBFile.writeFile(topology, state.getPositions(), fh)
        with open(nvt_chk, "wb") as fh:
            fh.write(sim.context.createCheckpoint())
        log.info("  NVT done → %s", cfg.omm_nvt_pdb.name)

    # ── NPT ───────────────────────────────────────────────────────────────────
    npt_chk = cfg.workdir / "npt.chk"
    if not _skip(cfg.omm_npt_pdb, "omm-npt"):
        npt_steps = cfg.steps_npt
        log.info("  NPT (%d steps = %d ps) …", npt_steps, int(cfg.ns_npt * 1000))
        # Must add barostat before creating the Simulation
        system_npt = openmm.XmlSerializer.deserialize(cfg.omm_system.read_text())
        system_npt.addForce(restraint)   # re-add restraint
        system_npt.addForce(openmm.MonteCarloBarostat(
            1.0 * unit.bar, cfg.temp * unit.kelvin, 25))
        integrator = openmm.LangevinMiddleIntegrator(
            cfg.temp * unit.kelvin, 1.0 / unit.picosecond,
            cfg.dt * unit.picoseconds)
        integrator.setRandomNumberSeed(43)
        sim = Simulation(topology, system_npt, integrator, _omm_platform(cfg))
        with open(nvt_chk, "rb") as fh:
            sim.context.loadCheckpoint(fh.read())
        sim.reporters.append(StateDataReporter(
            str(cfg.workdir / "npt_log.csv"), max(500, npt_steps // 100),
            step=True, time=True, temperature=True, potentialEnergy=True,
            volume=True, density=True))
        sim.step(npt_steps)
        state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
        with open(cfg.omm_npt_pdb, "w") as fh:
            PDBFile.writeFile(topology, state.getPositions(), fh)
        with open(npt_chk, "wb") as fh:
            fh.write(sim.context.createCheckpoint())
        log.info("  NPT done → %s", cfg.omm_npt_pdb.name)


# ─── OpenMM stage 5: production run ───────────────────────────────────────────

def omm_production(cfg: MDConfig) -> None:
    """
    OpenMM Stage 5 — Production MD.

    No position restraints. MonteCarloBarostat for NPT ensemble.
    Saves a DCD trajectory (md_traj.dcd) and energy CSV (omm_energy.csv).
    Checkpoints every ~1 ns so the run can be extended later.
    """
    _omm_check_import()
    log.info("[run/openmm] Production MD (%.1f ns, %d steps)",
             cfg.ns_prod, cfg.steps_prod)
    if _skip(cfg.omm_md_dcd, "omm-production"):
        return

    import openmm
    import openmm.unit as unit
    from openmm.app import (PDBFile, Simulation, DCDReporter,
                            StateDataReporter, CheckpointReporter)

    topology, _ = _omm_load_pdb(cfg.omm_topology)
    # Fresh system (no restraints) + barostat for NPT production
    system = openmm.XmlSerializer.deserialize(cfg.omm_system.read_text())
    system.addForce(openmm.MonteCarloBarostat(
        1.0 * unit.bar, cfg.temp * unit.kelvin, 25))

    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temp * unit.kelvin, 1.0 / unit.picosecond,
        cfg.dt * unit.picoseconds)
    integrator.setRandomNumberSeed(44)
    sim = Simulation(topology, system, integrator, _omm_platform(cfg))

    with open(cfg.workdir / "npt.chk", "rb") as fh:
        sim.context.loadCheckpoint(fh.read())

    sim.reporters.append(DCDReporter(str(cfg.omm_md_dcd), cfg.save_steps))
    sim.reporters.append(StateDataReporter(
        str(cfg.omm_energy_csv), cfg.save_steps,
        step=True, time=True, temperature=True,
        potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, volume=True, density=True,
        progress=True, totalSteps=cfg.steps_prod,
        remainingTime=True, speed=True))
    # Checkpoint every ~1 ns
    ckpt_interval = max(cfg.save_steps, round(1000 / cfg.dt))
    sim.reporters.append(CheckpointReporter(
        str(cfg.workdir / "md_running.chk"), ckpt_interval))

    log.info("  Running %.1f ns (dt=%.3f ps, save every %.0f ps) …",
             cfg.ns_prod, cfg.dt, cfg.save_ps)
    sim.step(cfg.steps_prod)

    state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
    with open(cfg.workdir / "md_final.pdb", "w") as fh:
        PDBFile.writeFile(topology, state.getPositions(), fh)
    log.info("  DCD → %s  |  Energy → %s",
             cfg.omm_md_dcd.name, cfg.omm_energy_csv.name)


# ─── OpenMM stage 6: analysis ─────────────────────────────────────────────────

def _omm_write_ref_and_fit(cfg: MDConfig, u, protein, ref_pdb: Path) -> None:
    """
    Write the frame-0 reference PDB (protein only) and a backbone-fitted
    protein-only DCD.  Both files use only protein atoms so MDAnalysis can
    load them together without an atom-count mismatch.
    """
    fit_dcd = cfg.omm_md_fit
    if ref_pdb.exists() and fit_dcd.exists():
        return

    import MDAnalysis as mda
    from MDAnalysis.analysis.align import alignto

    # Write frame-0 reference
    u.trajectory[0]
    if not ref_pdb.exists():
        with mda.Writer(str(ref_pdb)) as w:
            w.write(protein)
        log.info("  Reference PDB → %s", ref_pdb.name)

    # Write backbone-fitted protein-only DCD
    if not fit_dcd.exists():
        ref_u = mda.Universe(str(ref_pdb))
        n_frames = len(u.trajectory)
        log.info("  Fitting %d frames (backbone alignment) …", n_frames)
        with mda.Writer(str(fit_dcd), n_atoms=len(protein)) as writer:
            for _ts in u.trajectory:
                try:
                    alignto(u, ref_u, select="protein and backbone",
                            weights="mass")
                except Exception:
                    pass   # write unaligned frame if alignment fails
                writer.write(protein)
        log.info("  Fitted DCD → %s (%d protein atoms)", fit_dcd.name, len(protein))


def _parse_omm_energy_csv(csv_path: Path) -> Optional[np.ndarray]:
    """Parse the OpenMM StateDataReporter CSV into (n, 3): time_ps, epot, temp."""
    import csv as _csv
    rows = []
    try:
        with open(csv_path) as fh:
            reader = _csv.DictReader(fh)
            for row in reader:
                try:
                    t   = float(row.get("Time (ps)", 0) or 0)
                    ep  = float(row.get("Potential Energy (kJ/mole)", 0) or 0)
                    tem = float(row.get("Temperature (K)", 0) or 0)
                    rows.append([t, ep, tem])
                except (ValueError, TypeError):
                    continue
    except Exception as exc:
        log.warning("  Could not parse energy CSV: %s", exc)
        return None
    return np.array(rows) if rows else None


def omm_analyze(cfg: MDConfig) -> dict:
    """
    OpenMM Stage 6 — Trajectory analysis via MDAnalysis.

    Returns the same ``results`` dict shape as ``run_analysis()`` so that
    ``plot_results()``, ``write_pymol_script()``, and
    ``generate_html_report()`` all work unchanged.

    Unit convention (matches GROMACS XVG output that plot_results expects):
      rmsd[:, 0]  — time in ns          rmsd[:, 1]  — RMSD in nm
      rmsf[:, 0]  — residue number      rmsf[:, 1]  — RMSF in nm
      rg[:, 0]    — time in ps          rg[:, 1]    — Rg in nm
      energy[:, 0]— time in ps          energy[:, 1] — kJ/mol, energy[:, 2] — K
      hbonds[:, 0]— time in ns          hbonds[:, 1] — count
    """
    log.info("[analyze/openmm] Running MDAnalysis trajectory analyses")

    if not cfg.omm_md_dcd.exists():
        log.warning("  md_traj.dcd not found — skipping analysis")
        return {}

    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import rms as mda_rms
    except ImportError:
        log.warning("  MDAnalysis not installed — skipping (pip install MDAnalysis)")
        return {}

    try:
        u = mda.Universe(str(cfg.omm_topology), str(cfg.omm_md_dcd))
    except Exception as exc:
        log.warning("  Could not load trajectory: %s", exc)
        return {}

    protein  = u.select_atoms("protein")
    n_frames = len(u.trajectory)
    log.info("  Trajectory: %d frames, %d protein atoms", n_frames, len(protein))

    # Write reference PDB + fitted protein-only DCD (used by DSSP/PCA too)
    ref_pdb = cfg.workdir / "md_ref.pdb"
    _omm_write_ref_and_fit(cfg, u, protein, ref_pdb)

    # Re-load using fitted protein-only DCD for all analyses
    try:
        u_fit = mda.Universe(str(ref_pdb), str(cfg.omm_md_fit))
        prot  = u_fit.select_atoms("all")   # all atoms are protein in this universe
    except Exception as exc:
        log.warning("  Could not load fitted trajectory: %s — using raw", exc)
        u_fit = u
        prot  = protein

    results = {}

    # ── RMSD ──────────────────────────────────────────────────────────────────
    rmsd_npy = cfg.workdir / "omm_rmsd.npy"
    if not _skip(rmsd_npy, "omm-rmsd"):
        try:
            log.info("  RMSD …")
            R = mda_rms.RMSD(prot, select="backbone")
            R.run()
            # R.results.rmsd: (n_frames, 3) → [frame, time_ps, rmsd_Å]
            raw = R.results.rmsd
            arr = np.column_stack([raw[:, 1] / 1000.0,   # ps → ns
                                   raw[:, 2] / 10.0])    # Å → nm
            np.save(rmsd_npy, arr)
        except Exception as exc:
            log.warning("  RMSD failed: %s", exc)
            arr = np.empty((0, 2))
    else:
        arr = np.load(rmsd_npy)
    if len(arr):
        results["rmsd"] = arr
        log.info("  RMSD: %.3f ± %.3f Å", arr[:, 1].mean() * 10, arr[:, 1].std() * 10)

    # ── RMSF ──────────────────────────────────────────────────────────────────
    rmsf_npy = cfg.workdir / "omm_rmsf.npy"
    if not _skip(rmsf_npy, "omm-rmsf"):
        try:
            log.info("  RMSF …")
            ca = prot.select_atoms("name CA")
            rmsf_calc = mda_rms.RMSF(ca).run()
            arr_rmsf = np.column_stack([ca.resids,
                                        rmsf_calc.results.rmsf / 10.0])  # Å → nm
            np.save(rmsf_npy, arr_rmsf)
        except Exception as exc:
            log.warning("  RMSF failed: %s", exc)
            arr_rmsf = np.empty((0, 2))
    else:
        arr_rmsf = np.load(rmsf_npy)
    if len(arr_rmsf):
        results["rmsf"] = arr_rmsf

    # ── Radius of gyration ────────────────────────────────────────────────────
    rg_npy = cfg.workdir / "omm_rg.npy"
    if not _skip(rg_npy, "omm-rg"):
        try:
            log.info("  Radius of gyration …")
            times_ps, rg_nm = [], []
            for ts in u_fit.trajectory:
                times_ps.append(ts.time)
                rg_nm.append(prot.radius_of_gyration() / 10.0)  # Å → nm
            arr_rg = np.column_stack([times_ps, rg_nm])
            np.save(rg_npy, arr_rg)
        except Exception as exc:
            log.warning("  Rg failed: %s", exc)
            arr_rg = np.empty((0, 2))
    else:
        arr_rg = np.load(rg_npy)
    if len(arr_rg):
        results["rg"] = arr_rg

    # ── Energy from OpenMM CSV reporter ───────────────────────────────────────
    if cfg.omm_energy_csv.exists():
        log.info("  Parsing energy CSV …")
        energy_arr = _parse_omm_energy_csv(cfg.omm_energy_csv)
        if energy_arr is not None and len(energy_arr):
            results["energy"] = energy_arr

    # ── H-bonds ───────────────────────────────────────────────────────────────
    hbond_npy = cfg.workdir / "omm_hbonds.npy"
    if not _skip(hbond_npy, "omm-hbonds"):
        try:
            log.info("  H-bonds …")
            from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
            hba = HydrogenBondAnalysis(
                universe=u_fit,
                donors_sel="protein and name N",
                acceptors_sel="protein and name O* N*",
                d_a_cutoff=3.0,
                d_h_a_angle_cutoff=150.0,
            )
            hba.run()
            # count_by_time() → (n_frames, 2): [time_ps, count]
            counts = hba.count_by_time()
            arr_hb = np.column_stack([counts[:, 0] / 1000.0,  # ps → ns
                                      counts[:, 1]])
            np.save(hbond_npy, arr_hb)
        except Exception as exc:
            log.warning("  H-bond analysis failed (%s) — skipping", exc)
            arr_hb = np.empty((0, 2))
    else:
        arr_hb = np.load(hbond_npy)
    if len(arr_hb):
        results["hbonds"] = arr_hb

    # ── DSSP and PCA — reuse existing shared functions ─────────────────────────
    # md_ref.pdb (protein-only) + omm_md_fit (protein-only DCD) → atom counts match
    traj_for_mda = cfg.omm_md_fit if cfg.omm_md_fit.exists() else cfg.omm_md_dcd
    results["dssp"] = _compute_dssp(cfg, traj_for_mda)
    results["pca"]  = _compute_pca(cfg, traj_for_mda)

    return results


# ─── stage 6: trajectory processing + analysis ────────────────────────────────

def process_trajectory(cfg: MDConfig) -> None:
    """
    Pre-process the raw trajectory before analysis.

    Step 1 — Centre + PBC: protein is centred in the box and broken molecules
              (split across periodic boundaries) are made whole.
    Step 2 — Fit: rotational and translational motion of the protein backbone
              is removed so that structural fluctuations stand out clearly.

    GROMACS group indices (standard protein-in-water system):
      0=System  1=Protein  3=C-alpha  4=Backbone
    """
    log.info("[analyze] Processing trajectory (centre → fit)")

    if not cfg.md_tpr.exists():
        log.warning("  md.tpr not found — skipping trajectory processing")
        return

    if not _skip(cfg.md_center, "trjconv-center"):
        log.info("  Centering and repairing PBC …")
        gmx("trjconv",
            "-s", cfg.md_tpr, "-f", cfg.md_xtc, "-o", cfg.md_center,
            "-center", "-pbc", "mol", "-ur", "compact",
            stdin="1\n0\n",   # centre on Protein(1), write System(0)
            cwd=cfg.workdir)

    if not _skip(cfg.md_fit, "trjconv-fit"):
        log.info("  Fitting to backbone (remove rotation/translation) …")
        gmx("trjconv",
            "-s", cfg.md_tpr, "-f", cfg.md_center, "-o", cfg.md_fit,
            "-fit", "rot+trans",
            stdin="4\n1\n",   # fit on Backbone(4), write Protein(1)
            cwd=cfg.workdir)


def run_analysis(cfg: MDConfig) -> dict:
    """
    Stage 6 — Run all trajectory analyses.

    Returns a dict: {name: np.ndarray} for quantitative analyses.
    Each array has time in column 0 (units noted per tool below).

    Tools used:
      gmx rms    — Cα backbone RMSD vs initial structure (nm; col 0 in ns)
      gmx rmsf   — per-residue Cα RMSF over full trajectory (nm vs residue #)
      gmx gyrate — radius of gyration over time (nm; col 0 in ps)
      gmx energy — potential energy, temperature, pressure from EDR file
      gmx hbond  — intra-protein H-bond count over time
      MDAnalysis — DSSP secondary structure per residue per frame (optional)
    """
    log.info("[analyze] Running trajectory analyses")

    if not cfg.md_tpr.exists():
        log.warning("  md.tpr not found — skipping analysis")
        return {}

    traj = cfg.md_fit if cfg.md_fit.exists() else cfg.md_xtc
    d    = cfg.workdir
    results = {}

    # ── RMSD — backbone Cα vs starting structure ──────────────────────────────
    rmsd_xvg = d / "rmsd.xvg"
    data = _gmx_analysis(cfg, rmsd_xvg, "rmsd",
                         "rms", "-s", cfg.md_tpr, "-f", traj,
                         "-o", rmsd_xvg, "-tu", "ns",
                         stdin="4\n4\n")   # Backbone for fit + RMSD
    if data is not None and len(data):
        results["rmsd"] = data
        log.info("  RMSD: %.3f ± %.3f Å (mean ± SD over trajectory)",
                 data[:, 1].mean() * 10, data[:, 1].std() * 10)

    # ── RMSF — per-residue Cα fluctuations ───────────────────────────────────
    rmsf_xvg = d / "rmsf.xvg"
    data = _gmx_analysis(cfg, rmsf_xvg, "rmsf",
                         "rmsf", "-s", cfg.md_tpr, "-f", traj,
                         "-o", rmsf_xvg, "-res",
                         stdin="3\n")    # C-alpha
    if data is not None and len(data):
        results["rmsf"] = data

    # ── Radius of gyration ────────────────────────────────────────────────────
    gyrate_xvg = d / "gyrate.xvg"
    data = _gmx_analysis(cfg, gyrate_xvg, "gyrate",
                         "gyrate", "-s", cfg.md_tpr, "-f", traj,
                         "-o", gyrate_xvg,
                         stdin="1\n")    # Protein
    if data is not None and len(data):
        results["rg"] = data

    # ── Energy terms from EDR file ────────────────────────────────────────────
    if cfg.md_edr.exists():
        energy_xvg = d / "energy.xvg"
        data = _gmx_analysis(cfg, energy_xvg, "energy",
                             "energy", "-f", cfg.md_edr, "-o", energy_xvg,
                             stdin="Potential\nTemperature\nPressure\n\n")
        if data is not None and len(data):
            results["energy"] = data

    # ── Intra-protein hydrogen bond count ─────────────────────────────────────
    hbond_xvg = d / "hbond.xvg"
    data = _gmx_analysis(cfg, hbond_xvg, "hbond",
                         "hbond", "-s", cfg.md_tpr, "-f", traj,
                         "-num", hbond_xvg,
                         stdin="1\n1\n")  # Protein Protein
    if data is not None and len(data):
        results["hbonds"] = data

    # ── DSSP secondary structure (via MDAnalysis) ─────────────────────────────
    results["dssp"] = _compute_dssp(cfg, traj)

    # ── PCA of backbone dynamics (via MDAnalysis) ──────────────────────────────
    results["pca"] = _compute_pca(cfg, traj)

    return results


def _gmx_analysis(cfg: MDConfig, out: Path, label: str,
                  *gmx_args, stdin: str = "") -> Optional[np.ndarray]:
    """Run one GROMACS analysis tool and return parsed XVG data (or None)."""
    if not _skip(out, label):
        try:
            gmx(*gmx_args, stdin=stdin, cwd=cfg.workdir)
        except RuntimeError as e:
            log.warning("  [warn] %s — skipping: %s", label, e)
            return None
    if out.exists() and out.stat().st_size > 0:
        data, _ = parse_xvg(out)
        return data
    return None


def _compute_dssp(cfg: MDConfig, traj: Path) -> Optional[dict]:
    """Compute per-residue secondary structure using MDAnalysis DSSP."""
    dssp_npz = cfg.workdir / "dssp.npz"
    if dssp_npz.exists():
        npz = np.load(dssp_npz, allow_pickle=True)
        return dict(npz)

    try:
        import MDAnalysis as mda
        try:
            from MDAnalysis.analysis.dssp import DSSP
        except ImportError:
            try:
                from MDAnalysis.analysis.secondary_structure import DSSP
            except ImportError:
                log.warning("  DSSP not available in this MDAnalysis version — skipping")
                return None
    except ImportError:
        log.warning("  MDAnalysis not installed — skipping DSSP (pip install MDAnalysis)")
        return None

    if not traj.exists():
        return None

    # MDAnalysis may not support the newest GROMACS TPR format.
    # The fitted trajectory is protein-only, so we prefer md_ref.pdb as fallback topology.
    ref_pdb = cfg.workdir / "md_ref.pdb"
    if cfg.md_tpr.exists():
        top = cfg.md_tpr
    elif ref_pdb.exists():
        top = ref_pdb
    else:
        return None

    try:
        log.info("  DSSP via MDAnalysis …")
        try:
            u = mda.Universe(str(top), str(traj))
        except Exception:
            if top != ref_pdb and ref_pdb.exists():
                log.info("  TPR unsupported by MDAnalysis — falling back to md_ref.pdb")
                u = mda.Universe(str(ref_pdb), str(traj))
            else:
                raise
        # GROMACS names the C-terminal oxygens OC1/OC2 — rename to O/OXT so DSSP
        # sees equal counts of N, CA, C, O backbone atoms.
        for at in u.select_atoms("name OC1"):
            at.name = "O"
        for at in u.select_atoms("name OC2"):
            at.name = "OXT"
        protein = u.select_atoms("protein")
        ds      = DSSP(protein)
        ds.run()
        ss      = ds.results.dssp   # (n_frames, n_residues) char codes
        times   = np.array([ts.time for ts in u.trajectory]) / 1000.0  # → ns
        np.savez(dssp_npz, ss=ss, times=times)
        log.info("  DSSP: %d frames × %d residues", *ss.shape)
        return {"ss": ss, "times": times}
    except Exception as exc:
        log.warning("  DSSP failed: %s", exc)
        return None


def _compute_pca(cfg: MDConfig, traj: Path) -> Optional[dict]:
    """
    Compute Principal Component Analysis of backbone dynamics via MDAnalysis.

    Returns a dict with:
      variance    — fraction of variance per PC (first 20)
      cumvar      — cumulative variance (first 20)
      projected   — (n_frames, 3) projection of trajectory onto PC1-3
      times       — frame times in ns
      pc1_per_res — per-residue PC1 contribution magnitude (for B-factor colouring)
      pc1_resids  — residue numbers corresponding to pc1_per_res
    """
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis.pca import PCA
    except ImportError:
        log.warning("  MDAnalysis not installed — skipping PCA (pip install MDAnalysis)")
        return None

    pca_npz = cfg.workdir / "pca.npz"
    if pca_npz.exists():
        npz = np.load(pca_npz, allow_pickle=True)
        return {k: npz[k] for k in npz.files}

    if not traj.exists():
        return None

    # MDAnalysis may not support the newest GROMACS TPR format.
    # The fitted trajectory is protein-only, so we prefer md_ref.pdb as fallback topology.
    ref_pdb = cfg.workdir / "md_ref.pdb"
    if cfg.md_tpr.exists():
        top = cfg.md_tpr
    elif ref_pdb.exists():
        top = ref_pdb
    else:
        return None

    try:
        log.info("  PCA via MDAnalysis …")
        try:
            u = mda.Universe(str(top), str(traj))
        except Exception:
            if top != ref_pdb and ref_pdb.exists():
                log.info("  TPR unsupported by MDAnalysis — falling back to md_ref.pdb")
                u = mda.Universe(str(ref_pdb), str(traj))
            else:
                raise
        pc = PCA(u, select="backbone", align=True).run()

        # pc.results.variance contains raw eigenvalues (Å²); normalise to fractions
        raw_var   = pc.results.variance
        total_var = raw_var.sum() if raw_var.sum() > 0 else 1.0
        variance  = raw_var / total_var          # fraction per component
        cumvar    = np.cumsum(variance)

        # Project each frame onto first 3 PCs
        backbone  = u.select_atoms("backbone")
        projected = pc.transform(backbone, n_components=3)  # (n_frames, 3)

        # Frame times in ns
        times = np.array([ts.time for ts in u.trajectory]) / 1000.0

        # Per-residue PC1 contribution: RMS magnitude of x/y/z eigenvector components
        pc1_evec      = pc.results.p_components[0]           # (n_backbone_atoms * 3,)
        resids        = backbone.resids
        evec_3d       = pc1_evec.reshape(len(backbone), 3)
        pc1_per_atom  = np.sqrt((evec_3d ** 2).sum(axis=1))
        unique_res    = np.unique(resids)
        pc1_per_res   = np.array([pc1_per_atom[resids == r].mean() for r in unique_res])

        np.savez(pca_npz,
                 variance=variance[:20], cumvar=cumvar[:20],
                 projected=projected, times=times,
                 pc1_per_res=pc1_per_res, pc1_resids=unique_res)

        log.info("  PCA: PC1=%.1f%%  PC2=%.1f%%  (cumulative PC1+2=%.1f%%)",
                 variance[0]*100, variance[1]*100, cumvar[1]*100)

        return {"variance": variance[:20], "cumvar": cumvar[:20],
                "projected": projected, "times": times,
                "pc1_per_res": pc1_per_res, "pc1_resids": unique_res}
    except Exception as exc:
        log.warning("  PCA failed: %s", exc)
        return None


# ─── stage 7a: plots ──────────────────────────────────────────────────────────

def plot_results(cfg: MDConfig, results: dict) -> None:
    """
    Stage 7a — Save publication-quality matplotlib figures.

    Combined multi-panel summary figure + individual plot per analysis type:
      rmsd.png         Backbone RMSD vs time (Å)
      rmsf.png         Per-residue Cα RMSF bar chart (Å)
      rg.png           Radius of gyration vs time (Å)
      energy.png       Potential energy + temperature vs time
      hbonds.png       Intra-protein H-bond count vs time
      dssp.png         Secondary structure content (%) stacked area
      <stem>_summary.png  All panels combined in one figure
    """
    cfg.plots_dir.mkdir(exist_ok=True)
    log.info("[visualize] Plotting → %s/", cfg.plots_dir)

    PANELS = [
        ("rmsd",   "rmsd"   in results),
        ("rg",     "rg"     in results),
        ("rmsf",   "rmsf"   in results),
        ("hbonds", "hbonds" in results),
        ("energy", "energy" in results),
        ("dssp",   "dssp"   in results and results["dssp"] is not None),
        ("pca",    "pca"    in results and results["pca"] is not None),
    ]
    active = [name for name, ok in PANELS if ok]
    if not active:
        log.warning("  No analysis results — nothing to plot")
        return

    stem = cfg.pdb_in.stem

    # Combined figure
    nrows = len(active)
    fig   = plt.figure(figsize=(12, 3.5 * nrows))
    gs    = gridspec.GridSpec(nrows, 1, figure=fig, hspace=0.5)
    for i, name in enumerate(active):
        _render_panel(fig.add_subplot(gs[i]), name, results, cfg)
    fig.suptitle(f"MD Trajectory Analysis — {stem}", fontsize=13, y=1.01)
    combo = cfg.plots_dir / f"{stem}_summary.png"
    fig.savefig(combo, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("  Saved combined summary → %s", combo.name)

    # Individual panels
    for name in active:
        fig2, ax2 = plt.subplots(figsize=(10, 3.5))
        _render_panel(ax2, name, results, cfg)
        out = cfg.plots_dir / f"{stem}_{name}.png"
        fig2.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig2)
        log.info("  Saved  %s", out.name)

    # PCA scatter (separate multi-panel figure)
    if results.get("pca") is not None:
        _plot_pca_scatter(cfg, results["pca"])


def _render_panel(ax: plt.Axes, name: str, results: dict, cfg: MDConfig) -> None:
    """Draw one analysis panel onto ax."""
    stem = cfg.pdb_in.stem

    if name == "rmsd":
        d = results["rmsd"]
        ax.plot(d[:, 0], d[:, 1] * 10, color="#2077b4", lw=1)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Backbone RMSD (Å)")
        ax.set_title("Backbone RMSD vs. starting structure")
        _ax_style(ax)

    elif name == "rg":
        d = results["rg"]
        t = d[:, 0] / 1000   # ps → ns
        ax.plot(t, d[:, 1] * 10, color="#d62728", lw=1)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Rg (Å)")
        ax.set_title("Radius of gyration")
        _ax_style(ax)

    elif name == "rmsf":
        d = results["rmsf"]
        ax.bar(d[:, 0], d[:, 1] * 10, color="#2ca02c", width=0.8)
        ax.set_xlabel("Residue number")
        ax.set_ylabel("Cα RMSF (Å)")
        ax.set_title("Per-residue Cα flexibility")
        _ax_style(ax)

    elif name == "hbonds":
        d = results["hbonds"]
        t = d[:, 0]
        # hbond XVG time may be in ps — convert if range >> ns_prod
        if t.max() > cfg.ns_prod * 2:
            t = t / 1000
        ax.plot(t, d[:, 1], color="#9467bd", lw=1)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("H-bond count")
        ax.set_title("Intra-protein hydrogen bonds")
        _ax_style(ax)

    elif name == "energy":
        d = results["energy"]
        t = d[:, 0] / 1000   # ps → ns
        ax.plot(t, d[:, 1], color="#8c564b", lw=1, label="Potential (kJ/mol)")
        ax2 = ax.twinx()
        if d.shape[1] > 2:
            ax2.plot(t, d[:, 2], color="#e377c2", lw=0.8, ls="--",
                     label="Temperature (K)")
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Potential energy (kJ/mol)", color="#8c564b")
        ax2.set_ylabel("Temperature (K)", color="#e377c2")
        ax.set_title("Energy and temperature")
        ax.tick_params(axis="y", labelcolor="#8c564b")
        ax2.tick_params(axis="y", labelcolor="#e377c2")
        lines1, lbl1 = ax.get_legend_handles_labels()
        lines2, lbl2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, lbl1 + lbl2, fontsize=8, loc="upper right")

    elif name == "dssp":
        dd = results.get("dssp")
        if dd is not None:
            _render_dssp(ax, dd)
        else:
            ax.text(0.5, 0.5, "DSSP unavailable", ha="center", va="center",
                    transform=ax.transAxes)

    elif name == "pca":
        pd = results.get("pca")
        if pd is not None:
            _render_pca_scree(ax, pd)
        else:
            ax.text(0.5, 0.5, "PCA unavailable", ha="center", va="center",
                    transform=ax.transAxes)


def _render_dssp(ax: plt.Axes, dssp_data: dict) -> None:
    """Stacked area chart of secondary structure content over time."""
    ss    = dssp_data["ss"]          # (frames, residues) char array
    times = dssp_data["times"]       # ns
    helix  = np.mean(np.isin(ss, ["H", "G", "I"]), axis=1) * 100
    strand = np.mean(np.isin(ss, ["E", "B"]),       axis=1) * 100
    coil   = 100.0 - helix - strand
    ax.stackplot(times, helix, strand, coil,
                 labels=["α-helix / 3₁₀", "β-strand", "Coil / Turn"],
                 colors=["#d62728", "#1f77b4", "#aec7e8"], alpha=0.85)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Residue fraction (%)")
    ax.set_title("Secondary structure content")
    ax.set_ylim(0, 100)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", alpha=0.3)


def _render_pca_scree(ax: plt.Axes, pca_data: dict) -> None:
    """Bar chart of variance explained per PC + cumulative line."""
    var    = pca_data["variance"] * 100   # → percent
    cumvar = pca_data["cumvar"]   * 100
    n      = len(var)
    xs     = np.arange(1, n + 1)
    ax.bar(xs, var, color="#1f77b4", alpha=0.8, label="Per PC")
    ax2 = ax.twinx()
    ax2.plot(xs, cumvar, "o-", color="#ff7f0e", ms=4, lw=1.5,
             label="Cumulative")
    ax2.axhline(90, color="gray", ls="--", lw=0.8, alpha=0.7)
    ax2.set_ylim(0, 105)
    ax2.set_ylabel("Cumulative variance (%)", color="#ff7f0e")
    ax2.tick_params(axis="y", labelcolor="#ff7f0e")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained (%)", color="#1f77b4")
    ax.tick_params(axis="y", labelcolor="#1f77b4")
    ax.set_title("PCA scree plot — backbone dynamics")
    ax.set_xticks(xs)
    lines1, lbl1 = ax.get_legend_handles_labels()
    lines2, lbl2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, lbl1 + lbl2, fontsize=8, loc="center right")


def _plot_pca_scatter(cfg: MDConfig, pca_data: dict) -> None:
    """
    Two-panel PCA figure:
      Left:  PC1 vs PC2 scatter coloured by simulation time
      Right: PC1 and PC2 time series
    Saved as <stem>_pca_scatter.png
    """
    proj  = pca_data["projected"]   # (n_frames, 3) Å
    times = pca_data["times"]       # ns
    stem  = cfg.pdb_in.stem

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Left: scatter coloured by time
    sc = ax1.scatter(proj[:, 0], proj[:, 1], c=times,
                     cmap="viridis", s=6, alpha=0.7)
    plt.colorbar(sc, ax=ax1, label="Time (ns)")
    ax1.set_xlabel("PC1 (Å)")
    ax1.set_ylabel("PC2 (Å)")
    v0, v1 = pca_data["variance"][0]*100, pca_data["variance"][1]*100
    ax1.set_title(f"PC1 vs PC2  ({v0:.1f}% + {v1:.1f}% variance)")
    ax1.axhline(0, color="k", lw=0.5, alpha=0.3)
    ax1.axvline(0, color="k", lw=0.5, alpha=0.3)
    _ax_style(ax1)

    # Right: PC1 and PC2 vs time
    ax2.plot(times, proj[:, 0], color="#1f77b4", lw=0.8, label="PC1")
    ax2.plot(times, proj[:, 1], color="#ff7f0e", lw=0.8, label="PC2")
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel("Projection (Å)")
    ax2.set_title("PC1 & PC2 time series")
    ax2.legend(fontsize=8)
    _ax_style(ax2)

    fig.suptitle(f"PCA of backbone dynamics — {stem}", y=1.01)
    out = cfg.plots_dir / f"{stem}_pca_scatter.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("  Saved  %s", out.name)

    # Per-residue PC1 contribution bar chart
    if "pc1_per_res" in pca_data and "pc1_resids" in pca_data:
        fig3, ax3 = plt.subplots(figsize=(10, 3.5))
        resids = pca_data["pc1_resids"]
        contribs = pca_data["pc1_per_res"]
        ax3.bar(resids, contribs, color="#d62728", width=0.8, alpha=0.85)
        ax3.set_xlabel("Residue number")
        ax3.set_ylabel("PC1 eigenvector magnitude")
        ax3.set_title("Per-residue PC1 contribution (backbone dynamics driver)")
        _ax_style(ax3)
        out3 = cfg.plots_dir / f"{stem}_pca_residues.png"
        fig3.savefig(out3, dpi=150, bbox_inches="tight")
        plt.close(fig3)
        log.info("  Saved  %s", out3.name)


def _ax_style(ax: plt.Axes) -> None:
    ax.grid(axis="y", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ─── stage 7b: PyMOL script ───────────────────────────────────────────────────

def write_pymol_script(cfg: MDConfig, results: dict) -> None:
    """
    Stage 7b — Write a self-contained PyMOL session script (.pml).

    The script (open with: cd <workdir> && pymol <stem>_session.pml):
      - Loads the protein reference structure (first frame of fitted trajectory)
      - Loads the full fitted XTC trajectory onto that object
      - Applies secondary-structure-based colouring (helix=red, strand=yellow)
      - Optionally colours by RMSF (RMSF values written into B-factor column)
      - Sets up lighting, fog, and a smooth playback movie
      - Saves a PyMOL binary session (.pse) for later interactive use

    Requirements: PyMOL (open-source: https://github.com/schrodinger/pymol-open-source)
      macOS:   brew install pymol    or    conda install -c conda-forge pymol-open-source
      Ubuntu:  sudo apt install pymol
    """
    log.info("[visualize] Writing PyMOL script")

    # Extract frame 0 as reference PDB for structure loading
    ref_pdb = cfg.workdir / "md_ref.pdb"
    if not _skip(ref_pdb, "trjconv-ref") and cfg.md_tpr.exists():
        traj = cfg.md_fit if cfg.md_fit.exists() else cfg.md_xtc
        try:
            gmx("trjconv",
                "-s", cfg.md_tpr, "-f", traj, "-o", ref_pdb,
                "-dump", "0",
                stdin="1\n",   # output Protein
                cwd=cfg.workdir)
        except RuntimeError:
            log.warning("  Could not extract reference PDB")

    # Write RMSF and PC1 into B-factor columns for PyMOL colouring
    rmsf_pdb = _write_rmsf_pdb(cfg, results)
    pca_pdb  = _write_pc1_pdb(cfg, results)

    stem     = cfg.pdb_in.stem
    traj_rel = cfg.active_trajectory.name
    pml_path = cfg.workdir / f"{stem}_session.pml"

    lines = [
        f"# PyMOL session script — {stem}",
        f"# Generated by run_md.py",
        f"# Usage: cd {cfg.workdir.name}  &&  pymol {pml_path.name}",
        "",
        "# ── Display settings ─────────────────────────────────────────────",
        "bg_color white",
        "set ray_shadows, 0",
        "set antialias, 2",
        "set cartoon_fancy_helices, 1",
        "set cartoon_flat_sheets, 1",
        "set cartoon_loop_radius, 0.3",
        "set stick_radius, 0.15",
        "",
        "# ── Load structure + trajectory ───────────────────────────────────",
    ]

    if ref_pdb.exists():
        lines += [
            f"load {ref_pdb.name}, protein",
            f"load_traj {traj_rel}, protein, 1",
        ]
    else:
        lines += [
            f"# Reference PDB not found — load the GRO file manually:",
            f"# load {cfg.npt_gro.name}, protein",
            f"# load_traj {traj_rel}, protein, 1",
        ]

    lines += [
        "",
        "# ── Representations ──────────────────────────────────────────────",
        "hide everything, protein",
        "show cartoon, protein",
        "# show surface, protein",
        "# set transparency, 0.7",
        "",
        "# ── Colour by secondary structure (default view) ─────────────────",
        "color red,    protein and ss h",    # α-helix
        "color yellow, protein and ss s",    # β-strand
        "color cyan,   protein and ss l+''", # loop / coil
        "",
    ]

    if rmsf_pdb and rmsf_pdb.exists():
        lines += [
            "# ── RMSF-based colouring ─────────────────────────────────────",
            "# The file below has per-residue Cα RMSF in the B-factor column.",
            "# Uncomment these lines to switch from SS to RMSF colouring.",
            f"# load {rmsf_pdb.name}, rmsf_ref",
            "# alter protein, b=0.0",
            "# alter_state 1, rmsf_ref, b=b",
            "# cmd.alter('protein', 'b = cmd.get_model(\"rmsf_ref\", 1)"
            ".atom[int(index)-1].b')",
            "# spectrum b, blue_white_red, protein, minimum=0, maximum=3",
            "# delete rmsf_ref",
            "",
        ]

    if pca_pdb and pca_pdb.exists():
        lines += [
            "# ── PC1 eigenvector colouring (dominant motion driver) ────────",
            "# The file below has per-residue PC1 contribution in the B-factor column.",
            "# Uncomment to colour by PCA: white=rigid, red=highly mobile in PC1.",
            f"# load {pca_pdb.name}, pc1_ref",
            "# alter protein, b=0.0",
            "# cmd.alter('protein', 'b = cmd.get_model(\"pc1_ref\", 1)"
            ".atom[int(index)-1].b')",
            "# spectrum b, white_red, protein",
            "# delete pc1_ref",
            "",
        ]

    lines += [
        "# ── Camera ───────────────────────────────────────────────────────",
        "orient protein",
        "zoom protein, 5",
        "set depth_cue, 1",
        "set fog_start, 0.45",
        "",
        "# ── Trajectory playback ──────────────────────────────────────────",
        "# mset 1 -N plays states 1 through N sequentially (not 'repeat state 1 N times')",
        f"mset 1 -{cfg.n_frames}",
        "",
        "# ── Useful interactive commands ──────────────────────────────────",
        "# play              # play trajectory in viewer (press Escape to stop)",
        "# mplay             # play as movie",
        "# set movie_fps, 25 # playback speed",
        "# forward           # step one frame",
        "# png frame.png, dpi=300, ray=1  # save rendered image",
        "",
        "# ── Render movie (requires ffmpeg) ──────────────────────────────",
        "# import os; os.makedirs('frames', exist_ok=True)",
        "# mpng frames/frame",
        "# quit",
        "# # Then in shell:",
        f"# # ffmpeg -r 25 -i frames/frame%04d.png -c:v libx264 -crf 20 {stem}.mp4",
        "",
        "# ── Save PyMOL session for later ─────────────────────────────────",
        f"save {stem}.pse",
        "",
        "print 'Session saved. Type \"play\" to watch the trajectory.'",
    ]

    pml_path.write_text("\n".join(lines) + "\n")
    log.info("  PyMOL script → %s", pml_path.name)
    log.info("  Open with:   cd %s && pymol %s", cfg.workdir, pml_path.name)


def _write_rmsf_pdb(cfg: MDConfig, results: dict) -> Optional[Path]:
    """Write a copy of the reference PDB with RMSF (Å) in the B-factor column."""
    if "rmsf" not in results:
        return None
    ref_pdb = cfg.workdir / "md_ref.pdb"
    if not ref_pdb.exists():
        return None

    out = cfg.workdir / "rmsf_bfactor.pdb"
    if _skip(out, "rmsf_pdb"):
        return out

    rmsf_data = results["rmsf"]   # (n_res, 2): resnum, rmsf_nm
    rmsf_map  = {int(row[0]): row[1] * 10 for row in rmsf_data}  # nm → Å

    with open(ref_pdb) as fi, open(out, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try:
                    resi = int(line[22:26])
                    val  = min(rmsf_map.get(resi, 0.0), 99.99)
                    line = line[:60] + f"{val:6.2f}" + line[66:]
                except (ValueError, IndexError):
                    pass
            fo.write(line)

    log.info("  RMSF B-factor PDB → %s", out.name)
    return out


def _write_pc1_pdb(cfg: MDConfig, results: dict) -> Optional[Path]:
    """Write a copy of the reference PDB with PC1 eigenvector magnitude in the B-factor column."""
    pca_data = results.get("pca")
    if pca_data is None or "pc1_per_res" not in pca_data:
        return None
    ref_pdb = cfg.workdir / "md_ref.pdb"
    if not ref_pdb.exists():
        return None

    out = cfg.workdir / "pc1_bfactor.pdb"
    if _skip(out, "pc1_pdb"):
        return out

    resids   = pca_data["pc1_resids"]
    contribs = pca_data["pc1_per_res"]
    # Normalise to 0-99 range for B-factor column
    c_max    = contribs.max() if contribs.max() > 0 else 1.0
    pc1_map  = {int(r): float(c / c_max * 99) for r, c in zip(resids, contribs)}

    with open(ref_pdb) as fi, open(out, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try:
                    resi = int(line[22:26])
                    val  = min(pc1_map.get(resi, 0.0), 99.99)
                    line = line[:60] + f"{val:6.2f}" + line[66:]
                except (ValueError, IndexError):
                    pass
            fo.write(line)

    log.info("  PC1 B-factor PDB → %s", out.name)
    return out


# ─── HTML report ──────────────────────────────────────────────────────────────

def generate_html_report(cfg: MDConfig, results: dict) -> Path:
    """
    Generate a self-contained HTML report embedding all analysis plots as base64
    images plus a summary statistics table.

    Returns the path to the written HTML file.
    """
    import base64

    stem     = cfg.pdb_in.stem
    out_path = cfg.workdir / f"{stem}_report.html"
    log.info("[visualize] Writing HTML report → %s", out_path.name)

    def _img_tag(png: Path, alt: str = "", width: str = "100%") -> str:
        if not png.exists():
            return f'<p class="missing">Plot not available: {png.name}</p>'
        data = base64.b64encode(png.read_bytes()).decode()
        return (f'<img src="data:image/png;base64,{data}" '
                f'alt="{alt}" style="width:{width};max-width:900px;">')

    # ── Statistics cards ──────────────────────────────────────────────────────
    stats = []
    if "rmsd" in results:
        d = results["rmsd"][:, 1] * 10  # Å
        stats += [
            ("Backbone RMSD", f"{d.mean():.2f} ± {d.std():.2f} Å", "mean ± SD"),
            ("RMSD max", f"{d.max():.2f} Å", "peak deviation"),
        ]
    if "rg" in results:
        d = results["rg"][:, 1] * 10
        stats += [("Rg", f"{d.mean():.2f} ± {d.std():.2f} Å", "radius of gyration")]
    if "hbonds" in results:
        d = results["hbonds"][:, 1]
        stats += [("H-bonds", f"{d.mean():.1f} ± {d.std():.1f}", "intra-protein mean")]
    if "pca" in results and results["pca"] is not None:
        v = results["pca"]["variance"]
        stats += [
            ("PC1 variance", f"{v[0]*100:.1f}%", "dominant motion"),
            ("PC1+PC2", f"{(v[0]+v[1])*100:.1f}%", "cumulative"),
        ]

    cards_html = ""
    for label, value, note in stats:
        cards_html += f"""
        <div class="card">
          <div class="card-value">{value}</div>
          <div class="card-label">{label}</div>
          <div class="card-note">{note}</div>
        </div>"""

    # ── Per-analysis sections ─────────────────────────────────────────────────
    plots_dir = cfg.plots_dir

    sections = [
        ("Backbone RMSD",
         "Root-mean-square deviation of backbone Cα atoms vs the starting structure. "
         "Values below ~3 Å indicate the protein is stable.",
         plots_dir / f"{stem}_rmsd.png"),
        ("Radius of Gyration",
         "Overall protein compactness over time. A stable Rg indicates the protein "
         "does not unfold or collapse during the simulation.",
         plots_dir / f"{stem}_rg.png"),
        ("Per-residue RMSF",
         "Cα fluctuation per residue averaged over the trajectory. "
         "High RMSF residues are flexible; these are typically loops or termini.",
         plots_dir / f"{stem}_rmsf.png"),
        ("Hydrogen Bonds",
         "Number of intra-protein hydrogen bonds over time. "
         "A steady count indicates stable secondary structure.",
         plots_dir / f"{stem}_hbonds.png"),
        ("Energy & Temperature",
         "Potential energy and temperature over the production run. "
         "Both should be flat and stable after equilibration.",
         plots_dir / f"{stem}_energy.png"),
        ("Secondary Structure (DSSP)",
         "Fraction of residues in each secondary structure class over time. "
         "Requires MDAnalysis (pip install MDAnalysis).",
         plots_dir / f"{stem}_dssp.png"),
        ("PCA Scree Plot",
         "Variance explained by each principal component. "
         "The first few PCs capture the dominant conformational motions.",
         plots_dir / f"{stem}_pca.png"),
        ("PCA Phase Space",
         "Projection of the trajectory onto PC1 vs PC2 (coloured by time) "
         "and their time series. Clusters indicate distinct conformational states.",
         plots_dir / f"{stem}_pca_scatter.png"),
        ("PCA Per-residue Contribution",
         "Which residues drive the dominant motion (PC1). "
         "High bars mark the most dynamically active regions.",
         plots_dir / f"{stem}_pca_residues.png"),
    ]

    sections_html = ""
    for title, desc, png in sections:
        img = _img_tag(png, alt=title)
        sections_html += f"""
    <section>
      <h2>{title}</h2>
      <p>{desc}</p>
      {img}
    </section>"""

    # ── PyMOL instructions ────────────────────────────────────────────────────
    pml = cfg.workdir / f"{stem}_session.pml"
    pymol_block = ""
    if pml.exists():
        pymol_block = f"""
    <section>
      <h2>Interactive 3D Visualization</h2>
      <p>A PyMOL session script was generated. To open it:</p>
      <pre>cd {cfg.workdir}\npymol {pml.name}</pre>
      <p>Inside PyMOL, type <code>play</code> to watch the MD trajectory.
         Optional colouring by RMSF or PC1 magnitude is available — see comments
         inside the .pml file.</p>
    </section>"""

    # ── Simulation parameters table ───────────────────────────────────────────
    params = [
        ("PDB input",      cfg.pdb_in.name),
        ("Force field",    cfg.ff),
        ("Water model",    cfg.water),
        ("Temperature",    f"{cfg.temp} K"),
        ("Production",     f"{cfg.ns_prod} ns"),
        ("Frames",         str(cfg.n_frames)),
        ("Save interval",  f"{cfg.save_ps} ps"),
    ]
    rows = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in params)
    params_table = f"<table class='params'>{rows}</table>"

    # ── Assemble full HTML ────────────────────────────────────────────────────
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>MD Report — {stem}</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            margin: 0; background: #f5f5f5; color: #222; }}
    header {{ background: #1a3a5c; color: white; padding: 24px 40px; }}
    header h1 {{ margin: 0 0 4px; font-size: 1.6em; }}
    header p  {{ margin: 0; opacity: 0.8; font-size: 0.9em; }}
    main {{ max-width: 960px; margin: 32px auto; padding: 0 24px; }}
    .cards {{ display: flex; flex-wrap: wrap; gap: 16px; margin: 24px 0; }}
    .card {{ background: white; border-radius: 10px; padding: 16px 22px;
             box-shadow: 0 2px 8px rgba(0,0,0,.08); min-width: 140px; flex: 1; }}
    .card-value {{ font-size: 1.4em; font-weight: 700; color: #1a3a5c; }}
    .card-label {{ font-weight: 600; margin: 4px 0 2px; }}
    .card-note  {{ font-size: 0.8em; color: #666; }}
    section {{ background: white; border-radius: 10px; padding: 24px 28px;
               margin: 20px 0; box-shadow: 0 2px 8px rgba(0,0,0,.08); }}
    section h2 {{ margin: 0 0 10px; color: #1a3a5c; font-size: 1.15em; }}
    section p  {{ margin: 0 0 14px; color: #444; line-height: 1.5; }}
    img {{ border-radius: 6px; display: block; }}
    pre {{ background: #f0f0f0; padding: 12px 16px; border-radius: 6px;
           font-size: 0.85em; overflow-x: auto; }}
    code {{ background: #eee; padding: 2px 5px; border-radius: 3px; font-size: 0.9em; }}
    table.params {{ border-collapse: collapse; font-size: 0.9em; }}
    table.params th, table.params td {{ padding: 6px 14px; border-bottom: 1px solid #ddd;
                                        text-align: left; }}
    table.params th {{ width: 160px; color: #555; font-weight: 600; }}
    .missing {{ color: #999; font-style: italic; font-size: 0.85em; }}
    footer {{ text-align: center; color: #aaa; font-size: 0.8em; padding: 32px; }}
  </style>
</head>
<body>
<header>
  <h1>MD Trajectory Report — {stem}</h1>
  <p>Generated by run_md.py &nbsp;·&nbsp; {cfg.ns_prod} ns production MD
     &nbsp;·&nbsp; {cfg.ff} force field</p>
</header>
<main>

  <section>
    <h2>Simulation Parameters</h2>
    {params_table}
  </section>

  <section>
    <h2>Summary Statistics</h2>
    <div class="cards">{cards_html}
    </div>
  </section>

  {sections_html}

  {pymol_block}

</main>
<footer>Generated with run_md.py using {cfg.engine.upper()} &amp; MDAnalysis</footer>
</body>
</html>
"""

    out_path.write_text(html)
    log.info("  HTML report → %s", out_path.name)
    return out_path


# ─── CLI ──────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_md.py",
        description=textwrap.dedent("""\
            MD pipeline (GROMACS or OpenMM) — prepare → simulate → analyse → visualise.

            All stages run by default. Use --steps to run a subset:
              --steps prepare setup          build inputs only
              --steps minimize equil run     simulation only (inputs already built)
              --steps analyze visualize      post-process an existing run

            Engine selection:
              --engine gromacs   (default) requires GROMACS binary in PATH
              --engine openmm    requires: pip install openmm pdbfixer MDAnalysis
        """),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f"""\
            GROMACS force fields:
            {chr(10).join(f'  {k:20s}  {v}' for k, v in SUPPORTED_FORCEFIELDS.items())}

            OpenMM force fields:
            {chr(10).join(f'  {k:20s}  {v}' for k, v in OPENMM_FORCEFIELDS.items())}

            GROMACS water models:
            {chr(10).join(f'  {k:8s}  {v}' for k, v in WATER_MODELS.items())}

            OpenMM water models:
            {chr(10).join(f'  {k:8s}  {v}' for k, v in OPENMM_WATER_MODELS.items())}
        """),
    )
    p.add_argument("pdb", type=Path, help="Input PDB file")
    p.add_argument("--engine", choices=["gromacs", "openmm"], default=DEFAULT_ENGINE,
                   help="MD engine to use (default: gromacs)")

    p.add_argument("--steps", nargs="+", choices=ALL_STEPS, default=ALL_STEPS,
                   metavar="STEP",
                   help=f"Stages to run (default: all). Choices: {', '.join(ALL_STEPS)}")
    p.add_argument("--workdir", type=Path, default=None,
                   help="Output directory (default: <pdb_stem>_md/)")

    _all_ff    = list(SUPPORTED_FORCEFIELDS) + list(OPENMM_FORCEFIELDS)
    _all_water = list(WATER_MODELS) + list(OPENMM_WATER_MODELS)
    sim = p.add_argument_group("simulation parameters")
    sim.add_argument("--ff",    default=None,
                     help=(f"Force field. GROMACS default: {DEFAULT_FF}. "
                           f"OpenMM default: {DEFAULT_OMM_FF}. "
                           f"See epilog for all choices."))
    sim.add_argument("--water", default=None,
                     help=(f"Water model. GROMACS default: {DEFAULT_WATER}. "
                           f"OpenMM default: {DEFAULT_OMM_WATER}. "
                           f"See epilog for all choices."))
    sim.add_argument("--ns",    type=float, default=DEFAULT_NS,
                     help=f"Production MD length in nanoseconds (default: {DEFAULT_NS})")
    sim.add_argument("--temp",  type=float, default=DEFAULT_TEMP,
                     help=f"Temperature in K (default: {DEFAULT_TEMP})")
    sim.add_argument("--conc",  type=float, default=DEFAULT_CONC,
                     help=f"NaCl concentration in mol/L (default: {DEFAULT_CONC})")
    sim.add_argument("--box",   type=float, default=DEFAULT_BOX,
                     help=f"Min protein-to-box-edge in nm (default: {DEFAULT_BOX})")
    sim.add_argument("--save-ps", type=float, default=DEFAULT_SAVEPS,
                     help=f"Trajectory save interval in ps (default: {DEFAULT_SAVEPS})")
    sim.add_argument("--nvt-ps", type=float, default=100.0,
                     help="NVT equilibration length in ps (default: 100)")
    sim.add_argument("--npt-ps", type=float, default=100.0,
                     help="NPT equilibration length in ps (default: 100)")

    hw = p.add_argument_group("hardware")
    hw.add_argument("--ncores", type=int, default=0,
                    help="CPU threads for mdrun (default: auto-detect)")
    hw.add_argument("--gpu",    action="store_true",
                    help="Enable GPU offloading in mdrun")

    p.add_argument("-v", "--verbose", action="store_true",
                   help="Show full GROMACS command output")
    return p


def main():
    global GMX

    args = _build_parser().parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    engine = args.engine

    # Resolve per-engine defaults for --ff and --water
    if args.ff is None:
        args.ff = DEFAULT_FF if engine == "gromacs" else DEFAULT_OMM_FF
    if args.water is None:
        args.water = DEFAULT_WATER if engine == "gromacs" else DEFAULT_OMM_WATER

    # Validate force field / water model choices for the chosen engine
    if engine == "gromacs":
        if args.ff not in SUPPORTED_FORCEFIELDS:
            sys.exit(f"[error] Unknown GROMACS force field: {args.ff}\n"
                     f"Choices: {', '.join(SUPPORTED_FORCEFIELDS)}")
        if args.water not in WATER_MODELS:
            sys.exit(f"[error] Unknown GROMACS water model: {args.water}\n"
                     f"Choices: {', '.join(WATER_MODELS)}")
        GMX = _find_gmx()
        log.info("GROMACS: %s", shutil.which(GMX))
    else:
        _omm_check_import()
        import openmm
        log.info("OpenMM:  %s", openmm.__version__)
        if args.ff not in OPENMM_FORCEFIELDS:
            sys.exit(f"[error] Unknown OpenMM force field: {args.ff}\n"
                     f"Choices: {', '.join(OPENMM_FORCEFIELDS)}")
        if args.water not in OPENMM_WATER_MODELS:
            sys.exit(f"[error] Unknown OpenMM water model: {args.water}\n"
                     f"Choices: {', '.join(OPENMM_WATER_MODELS)}")

    pdb_in = args.pdb.resolve()
    if not pdb_in.exists():
        sys.exit(f"[error] PDB file not found: {pdb_in}")

    workdir = args.workdir or pdb_in.parent / f"{pdb_in.stem}_md"
    workdir.mkdir(parents=True, exist_ok=True)

    cfg = MDConfig(
        pdb_in   = pdb_in,
        workdir  = workdir,
        ff       = args.ff,
        water    = args.water,
        ns_prod  = args.ns,
        temp     = args.temp,
        ion_conc = args.conc,
        box_dist = args.box,
        save_ps  = args.save_ps,
        ns_nvt   = args.nvt_ps / 1000,
        ns_npt   = args.npt_ps / 1000,
        ncores   = args.ncores,
        gpu      = args.gpu,
        steps    = args.steps,
        engine   = engine,
    )

    log.info("═" * 62)
    log.info("  MD Pipeline  —  %s", pdb_in.name)
    log.info("  Engine:       %s", engine.upper())
    log.info("  Working dir:  %s", workdir)
    log.info("  Force field:  %s   Water: %s", cfg.ff, cfg.water)
    log.info("  Production:   %.0f ns at %.0f K  (dt = %.3f ps)",
             cfg.ns_prod, cfg.temp, cfg.dt)
    log.info("  Steps:        %s", ", ".join(cfg.steps))
    log.info("═" * 62)

    try:
        if "prepare" in cfg.steps:
            prepare_structure(cfg)

        if engine == "gromacs":
            if "setup"    in cfg.steps: setup_topology(cfg)
            if "minimize" in cfg.steps: energy_minimize(cfg)
            if "equil"    in cfg.steps: equilibrate(cfg)
            if "run"      in cfg.steps: production_run(cfg)
        else:  # openmm
            if "setup"    in cfg.steps: omm_setup(cfg)
            if "minimize" in cfg.steps: omm_minimize(cfg)
            if "equil"    in cfg.steps: omm_equilibrate(cfg)
            if "run"      in cfg.steps: omm_production(cfg)

        results = {}
        if "analyze" in cfg.steps or "visualize" in cfg.steps:
            if engine == "gromacs":
                process_trajectory(cfg)
                results = run_analysis(cfg)
            else:
                results = omm_analyze(cfg)

        if "visualize" in cfg.steps:
            plot_results(cfg, results)
            write_pymol_script(cfg, results)
            generate_html_report(cfg, results)

    except RuntimeError as exc:
        log.error("\n[PIPELINE FAILED] %s", exc)
        log.error("Logs and intermediate files are in: %s/", workdir)
        log.error("Fix the problem and re-run — completed steps will be skipped.")
        sys.exit(1)
    except KeyboardInterrupt:
        log.warning("\n[Interrupted] Re-run the same command to resume from this point.")
        sys.exit(130)

    log.info("")
    log.info("═" * 62)
    log.info("  Pipeline complete!")
    log.info("  Outputs in: %s/", workdir)
    log.info("  Plots in:   %s/plots/", workdir)
    pml = workdir / f"{pdb_in.stem}_session.pml"
    if pml.exists():
        log.info("  PyMOL:  cd %s && pymol %s", workdir, pml.name)
    html = workdir / f"{pdb_in.stem}_report.html"
    if html.exists():
        log.info("  Report: open %s", html)
    log.info("═" * 62)


if __name__ == "__main__":
    main()
