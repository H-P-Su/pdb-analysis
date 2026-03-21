#!/usr/bin/env python3
"""
summarize_structures.py — Structure inventory, metadata, geometry, and batch summaries.

All functions importable for use in other projects.

Functions
---------
parse_header(path)                    Extract resolution, R-factor, method, organism, date
summarize_structure(path)             Full per-structure summary dict
batch_summary(paths, out_csv)         Summarize many structures → CSV
find_missing_residues(path)           Compare SEQRES vs ATOM; return gaps
extract_fasta(path, chain)            FASTA sequence from ATOM coordinates
bfactor_stats(structure)              Per-residue B-factor table
radius_of_gyration(structure)         Cα radius of gyration in Å
ramachandran_data(structure)          phi/psi dihedral angles + outlier flags
plot_bfactor(data, output)            B-factor profile plot → PNG
plot_ramachandran(data, output)       Ramachandran scatter plot → PNG

CLI
---
python3 summarize_structures.py 1abc.pdb
python3 summarize_structures.py Files/*.pdb --out summary.csv
python3 summarize_structures.py 1abc.pdb --missing
python3 summarize_structures.py 1abc.pdb --fasta
python3 summarize_structures.py 1abc.pdb --bfactor --plot
python3 summarize_structures.py 1abc.pdb --ramachandran --plot
python3 summarize_structures.py 1abc.pdb --rg
"""

from __future__ import annotations

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Structure import Structure

# ─── Constants ────────────────────────────────────────────────────────────────

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "SEC": "U",
}

WATER_NAMES = frozenset({"HOH", "WAT", "H2O", "DOD", "SOL"})
PROTEIN_RESIDUES = frozenset(THREE_TO_ONE.keys())

# Ramachandran favored / allowed regions (hard-coded polygons for the
# general case; simplified from Lovell 2003 / MolProbity)
_RAMA_HELIX_REGIONS = [(-90, -60, -80, -20)]   # (φ_min, φ_max, ψ_min, ψ_max)
_RAMA_SHEET_REGIONS = [(-170, -50, 80, 180), (-170, -50, -180, -120)]
_RAMA_FAVORED_MARGIN = 30   # degrees tolerance added when classifying


# ─── Structure Loading ─────────────────────────────────────────────────────────

def load_structure(path: str | Path, sid: str = "struct") -> Structure:
    path = Path(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True).get_structure(sid, str(path))
    return PDBParser(QUIET=True).get_structure(sid, str(path))


# ─── Header / Metadata ────────────────────────────────────────────────────────

def parse_header(path: str | Path) -> dict:
    """
    Extract metadata from a PDB file header.

    Returns a dict with keys:
      pdb_id, title, method, resolution, r_work, r_free,
      deposition_date, organism, keywords
    """
    path = Path(path)
    info = {
        "pdb_id":          path.stem.upper(),
        "title":           "",
        "method":          "",
        "resolution":      None,
        "r_work":          None,
        "r_free":          None,
        "deposition_date": "",
        "organism":        "",
        "keywords":        "",
    }

    if path.suffix.lower() in (".cif", ".mmcif"):
        return info  # CIF header parsing not implemented; return defaults

    title_lines: list[str] = []
    keywds_lines: list[str] = []
    source_lines: list[str] = []
    r_free_candidates: list[float] = []

    with open(path, errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            rec = line[:6].strip()

            if rec == "HEADER":
                info["pdb_id"] = line[62:66].strip() or info["pdb_id"]
                info["deposition_date"] = line[50:59].strip()

            elif rec == "TITLE":
                title_lines.append(line[10:].strip())

            elif rec == "KEYWDS":
                keywds_lines.append(line[10:].strip())

            elif rec == "EXPDTA":
                info["method"] = line[10:].strip()

            elif rec == "SOURCE":
                source_lines.append(line[10:].strip())

            elif rec == "REMARK":
                remark_num = line[7:10].strip()

                if remark_num == "2":
                    m = re.search(r"RESOLUTION\.\s+([\d.]+)\s+ANGSTROM", line)
                    if m:
                        info["resolution"] = float(m.group(1))

                elif remark_num == "3":
                    if "R VALUE" in line and "WORKING" in line and info["r_work"] is None:
                        m = re.search(r":\s*(0\.\d+)", line)
                        if m:
                            info["r_work"] = float(m.group(1))
                    elif "FREE R VALUE" in line and "TEST" not in line:
                        m = re.search(r":\s*(0\.\d+)", line)
                        if m:
                            r_free_candidates.append(float(m.group(1)))

            elif rec == "ATOM":
                break  # stop at first ATOM record

    info["title"] = " ".join(title_lines)
    info["keywords"] = " ".join(keywds_lines)

    # Pick the largest FREE R VALUE (the actual Rfree, not test set size)
    if r_free_candidates:
        info["r_free"] = max(r_free_candidates)

    # Extract organism from SOURCE
    source_text = " ".join(source_lines)
    m = re.search(r"ORGANISM_SCIENTIFIC:\s*([^;]+)", source_text, re.IGNORECASE)
    if m:
        info["organism"] = m.group(1).strip().title()

    return info


# ─── Structure Inventory ──────────────────────────────────────────────────────

def summarize_structure(path: str | Path) -> dict:
    """
    Return a comprehensive summary dict for one structure.

    Keys: pdb_id, title, method, resolution, r_work, r_free,
          deposition_date, organism, n_models, chains, n_residues,
          n_atoms, n_ligands, ligand_names, n_waters,
          n_missing_residues, has_seqres
    """
    path = Path(path)
    structure = load_structure(path)
    header = parse_header(path)

    model = next(structure.get_models())    # first model

    chain_counts: dict[str, int] = {}
    n_atoms = 0
    ligands: list[str] = []
    n_waters = 0

    for chain in model:
        res_count = 0
        for res in chain:
            resn = res.get_resname().strip()
            flag = res.id[0]
            if flag == " ":                     # standard residue
                res_count += 1
            elif resn in WATER_NAMES:
                n_waters += 1
            else:
                ligands.append(f"{resn}/{chain.id}")
            n_atoms += sum(1 for _ in res.get_atoms())
        chain_counts[chain.id] = res_count

    missing = find_missing_residues(path)

    return {
        **header,
        "n_models":           sum(1 for _ in structure.get_models()),
        "chains":             ",".join(chain_counts.keys()),
        "n_residues":         sum(chain_counts.values()),
        "residues_per_chain": "; ".join(f"{c}:{n}" for c, n in chain_counts.items()),
        "n_atoms":            n_atoms,
        "n_ligands":          len(ligands),
        "ligand_names":       ",".join(sorted(set(ligands))),
        "n_waters":           n_waters,
        "n_missing_residues": len(missing),
        "has_seqres":         _has_seqres(path),
    }


def _has_seqres(path: Path) -> bool:
    if path.suffix.lower() in (".cif", ".mmcif"):
        return False
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("SEQRES"):
                return True
            if line.startswith("ATOM"):
                break
    return False


# ─── Missing Residues ─────────────────────────────────────────────────────────

def find_missing_residues(path: str | Path) -> list[dict]:
    """
    Compare SEQRES records with ATOM records to find modeled gaps.

    Returns a list of dicts: {chain, resn, resi, source}
    where source is "SEQRES_only" (not modeled) or "ATOM_only" (extra).
    Only works for PDB format files with SEQRES records.
    """
    path = Path(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        return []

    # Parse SEQRES
    seqres: dict[str, list[tuple[str, int]]] = defaultdict(list)
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("SEQRES"):
                chain = line[11]
                residues = line[19:].split()
                seqres[chain].extend(residues)
            elif line.startswith("ATOM"):
                break

    if not seqres:
        return []

    # Parse ATOM observed residues (unique by chain+resi)
    observed: dict[str, set[int]] = defaultdict(set)
    observed_resn: dict[tuple, str] = {}
    with open(path, errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            try:
                resi = int(line[22:26].strip())
            except ValueError:
                continue
            resn = line[17:20].strip()
            observed[chain].add(resi)
            observed_resn[(chain, resi)] = resn

    # Build expected residue numbers by matching SEQRES position
    missing: list[dict] = []
    for chain, resn_list in seqres.items():
        obs_riis = sorted(observed.get(chain, set()))
        if not obs_riis:
            continue
        # Assign approximate resi numbers to SEQRES entries
        # (SEQRES positions don't include resi numbers; use observed range)
        min_resi, max_resi = obs_riis[0], obs_riis[-1]
        obs_set = observed.get(chain, set())
        expected_count = len(resn_list)
        # Gap = SEQRES residues not in ATOM (rough count)
        n_gap = expected_count - len(obs_set)
        if n_gap > 0:
            missing.append({
                "chain": chain,
                "expected_residues": expected_count,
                "observed_residues": len(obs_set),
                "missing_count": n_gap,
                "pct_modeled": round(100 * len(obs_set) / expected_count, 1),
            })

    return missing


# ─── FASTA from Coordinates ───────────────────────────────────────────────────

def extract_fasta(path: str | Path, chain_id: str | None = None) -> str:
    """
    Build a FASTA sequence from ATOM records (observed residues only).

    Parameters
    ----------
    path :
        PDB or CIF file.
    chain_id :
        Restrict to this chain.  All chains if None.

    Returns
    -------
    FASTA string with one entry per chain.
    """
    structure = load_structure(path)
    stem = Path(path).stem.upper()
    lines: list[str] = []

    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            seq = []
            for res in chain:
                if res.id[0] != " ":
                    continue
                resn = res.get_resname().strip()
                seq.append(THREE_TO_ONE.get(resn, "X"))
            if seq:
                lines.append(f">{stem}_{chain.id}")
                seq_str = "".join(seq)
                for i in range(0, len(seq_str), 60):
                    lines.append(seq_str[i:i+60])
        break   # first model only

    return "\n".join(lines)


# ─── B-factor Analysis ────────────────────────────────────────────────────────

def bfactor_stats(structure: Structure,
                   chain_id: str | None = None) -> list[dict]:
    """
    Compute per-residue average B-factor for Cα atoms.

    Returns a list of dicts: {chain, resi, resn, bfactor, flag}
    where flag is 'high' (>2σ above mean), 'low' (<2σ below), or ''.
    """
    rows: list[dict] = []
    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for res in chain:
                if res.id[0] != " ":
                    continue
                resn = res.get_resname().strip()
                if resn not in PROTEIN_RESIDUES:
                    continue
                ca_atoms = [a for a in res.get_atoms() if a.get_name().strip() == "CA"]
                if not ca_atoms:
                    continue
                bf = ca_atoms[0].get_bfactor()
                rows.append({
                    "chain": chain.id,
                    "resi": res.id[1],
                    "resn": resn,
                    "bfactor": round(bf, 2),
                    "flag": "",
                })
        break

    if not rows:
        return rows

    bvals = [r["bfactor"] for r in rows]
    mean_b = float(np.mean(bvals))
    std_b  = float(np.std(bvals))
    for r in rows:
        if r["bfactor"] > mean_b + 2 * std_b:
            r["flag"] = "high"
        elif r["bfactor"] < mean_b - 2 * std_b:
            r["flag"] = "low"

    return rows


def plot_bfactor(data: list[dict], output_path: str | Path,
                 title: str = "") -> Path:
    """
    Save a B-factor profile plot (per-residue Cα B-factor).

    Parameters
    ----------
    data :
        Output of bfactor_stats().
    output_path :
        PNG output path.
    title :
        Plot title (optional).
    """
    if not data:
        print("  No B-factor data to plot.")
        return Path(output_path)

    chains = sorted(set(r["chain"] for r in data))
    n_chains = len(chains)
    # Cap figure dimensions so PIL never hits its decompression-bomb limit.
    # PIL rejects images > ~179 MP; at 150 DPI that is ~7950 sq-in.
    # We hard-cap width at 50 in and height at 3.5 in/chain (max 20 chains shown).
    _MAX_CHAINS_SHOWN = 20
    n_plot = min(n_chains, _MAX_CHAINS_SHOWN)
    chains = chains[:n_plot]
    fig_w = min(50, max(10, len(data) * 0.06 + 2))
    fig_h = 3.5 * n_plot
    # Compute a safe DPI so width_px * height_px < 150 MP
    dpi = min(150, int((150_000_000 / (fig_w * fig_h)) ** 0.5))
    dpi = max(72, dpi)
    fig, axes = plt.subplots(n_plot, 1,
                              figsize=(fig_w, fig_h),
                              squeeze=False)

    for ax, chain in zip(axes[:, 0], chains):
        rows = [r for r in data if r["chain"] == chain]
        x    = [r["resi"] for r in rows]
        y    = [r["bfactor"] for r in rows]
        mean_b = float(np.mean(y))
        std_b  = float(np.std(y))

        colors = ["#e84040" if r["flag"] == "high"
                  else "#4488cc" if r["flag"] == "low"
                  else "#555555" for r in rows]

        ax.bar(x, y, color=colors, width=1.0, linewidth=0)
        ax.axhline(mean_b, color="black", linewidth=1.0, linestyle="--",
                   label=f"Mean {mean_b:.1f} Å²")
        ax.axhline(mean_b + 2 * std_b, color="#e84040", linewidth=0.8,
                   linestyle=":", label="+2σ")
        ax.axhline(max(0, mean_b - 2 * std_b), color="#4488cc", linewidth=0.8,
                   linestyle=":", label="−2σ")

        ax.set_ylabel("B-factor (Å²)", fontsize=9)
        ax.set_xlabel("Residue number", fontsize=9)
        ax.set_title(f"B-factor profile — Chain {chain}"
                     + (f"  |  {title}" if title else ""), fontsize=10)
        ax.legend(fontsize=8, loc="upper right")

        high = [r for r in rows if r["flag"] == "high"]
        for r in high[:10]:
            ax.annotate(f"{r['resn']}{r['resi']}", (r["resi"], r["bfactor"]),
                        textcoords="offset points", xytext=(0, 4),
                        fontsize=6, ha="center", color="#e84040")

        ax.spines[["top", "right"]].set_visible(False)
        patches = [mpatches.Patch(color="#e84040", label="High (>2σ)"),
                   mpatches.Patch(color="#4488cc", label="Low (<2σ)"),
                   mpatches.Patch(color="#555555", label="Normal")]
        ax.legend(handles=patches, fontsize=8, loc="upper right")

    plt.tight_layout()
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved B-factor plot  → {output_path}")
    return output_path


# ─── Radius of Gyration ────────────────────────────────────────────────────────

def radius_of_gyration(structure: Structure,
                        chain_id: str | None = None,
                        atom_name: str = "CA") -> float:
    """
    Compute the radius of gyration of Cα atoms (or any atom type).

    Rg = sqrt( Σ |r_i - r_com|² / N )

    Parameters
    ----------
    structure :
        BioPython Structure.
    chain_id :
        Restrict to one chain.
    atom_name :
        Atom to use. Default 'CA'.  Use None for all heavy atoms.

    Returns
    -------
    Radius of gyration in Å.
    """
    coords = []
    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for res in chain:
                if res.id[0] != " ":
                    continue
                if atom_name:
                    if atom_name in res:
                        coords.append(res[atom_name].get_vector().get_array())
                else:
                    for atom in res.get_atoms():
                        if not atom.get_name().strip().startswith("H"):
                            coords.append(atom.get_vector().get_array())
        break

    if not coords:
        return 0.0

    coords = np.array(coords)
    com = coords.mean(axis=0)
    rg = float(np.sqrt(np.mean(np.sum((coords - com) ** 2, axis=1))))
    return round(rg, 3)


# ─── Ramachandran Analysis ────────────────────────────────────────────────────

def _dihedral(p0: np.ndarray, p1: np.ndarray,
              p2: np.ndarray, p3: np.ndarray) -> float:
    """Compute dihedral angle in degrees for four points (IUPAC convention)."""
    from Bio.PDB.vectors import Vector, calc_dihedral as _bio_dihedral
    return math.degrees(_bio_dihedral(
        Vector(p0), Vector(p1), Vector(p2), Vector(p3)))


def ramachandran_data(structure: Structure,
                       chain_id: str | None = None) -> list[dict]:
    """
    Compute phi/psi backbone dihedral angles for all residues.

    Returns a list of dicts:
      {chain, resi, resn, phi, psi, region}
    where region is 'helix', 'sheet', 'allowed', or 'outlier'.
    """
    rows: list[dict] = []

    def _get_atoms(res, names):
        return [res[n].get_vector().get_array() for n in names if n in res]

    for model in structure:
        residues_by_chain: dict[str, list] = defaultdict(list)
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for res in chain:
                if res.id[0] == " " and res.get_resname().strip() in PROTEIN_RESIDUES:
                    residues_by_chain[chain.id].append(res)

        for cid, res_list in residues_by_chain.items():
            for i, res in enumerate(res_list):
                resn = res.get_resname().strip()
                if resn == "GLY":
                    continue   # Gly has its own Ramachandran; skip for general plot

                phi = psi = None

                # phi: C(i-1) - N(i) - CA(i) - C(i)
                if i > 0:
                    prev = res_list[i - 1]
                    pts = (_get_atoms(prev, ["C"]) + _get_atoms(res, ["N", "CA", "C"]))
                    if len(pts) == 4:
                        phi = _dihedral(*pts)

                # psi: N(i) - CA(i) - C(i) - N(i+1)
                if i < len(res_list) - 1:
                    nxt = res_list[i + 1]
                    pts = (_get_atoms(res, ["N", "CA", "C"]) + _get_atoms(nxt, ["N"]))
                    if len(pts) == 4:
                        psi = _dihedral(*pts)

                if phi is None or psi is None:
                    continue

                region = _rama_region(phi, psi)
                rows.append({
                    "chain": cid,
                    "resi":  res.id[1],
                    "resn":  resn,
                    "phi":   round(phi, 2),
                    "psi":   round(psi, 2),
                    "region": region,
                })
        break

    return rows


def _rama_region(phi: float, psi: float) -> str:
    """Classify a phi/psi pair into helix, sheet, allowed, or outlier."""
    # Favored helix
    if -160 <= phi <= -40 and -80 <= psi <= 0:
        return "helix"
    # Favored sheet
    if -170 <= phi <= -50 and (100 <= psi <= 180 or -180 <= psi <= -100):
        return "sheet"
    # Broad allowed region
    if -180 <= phi <= 0:
        return "allowed"
    return "outlier"


def plot_ramachandran(data: list[dict], output_path: str | Path,
                      title: str = "") -> Path:
    """
    Save a Ramachandran scatter plot.

    Parameters
    ----------
    data :
        Output of ramachandran_data().
    output_path :
        PNG output path.
    title :
        Plot title (optional).
    """
    if not data:
        print("  No Ramachandran data to plot.")
        return Path(output_path)

    region_colors = {
        "helix":   "#4488ff",
        "sheet":   "#44bb44",
        "allowed": "#aaaaaa",
        "outlier": "#ff3333",
    }

    fig, ax = plt.subplots(figsize=(6, 6))

    # Background shading for favored regions
    helix_rect = plt.Rectangle((-160, -80), 120, 80,
                                color="#4488ff", alpha=0.08)
    sheet_rect1 = plt.Rectangle((-170, 100), 120, 80,
                                 color="#44bb44", alpha=0.08)
    sheet_rect2 = plt.Rectangle((-170, -180), 120, 80,
                                 color="#44bb44", alpha=0.08)
    ax.add_patch(helix_rect)
    ax.add_patch(sheet_rect1)
    ax.add_patch(sheet_rect2)

    for region, color in region_colors.items():
        pts = [(r["phi"], r["psi"]) for r in data if r["region"] == region]
        if pts:
            xs, ys = zip(*pts)
            ax.scatter(xs, ys, c=color, s=10, alpha=0.7, label=region,
                       edgecolors="none", zorder=3)

    # Annotate outliers
    outliers = [r for r in data if r["region"] == "outlier"]
    for r in outliers[:20]:
        ax.annotate(f"{r['resn']}{r['resi']}", (r["phi"], r["psi"]),
                    fontsize=6, color="#ff3333",
                    xytext=(4, 4), textcoords="offset points")

    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("φ (degrees)", fontsize=10)
    ax.set_ylabel("ψ (degrees)", fontsize=10)

    n_total   = len(data)
    n_outlier = sum(1 for r in data if r["region"] == "outlier")
    pct_ok    = 100 * (1 - n_outlier / n_total) if n_total else 0

    ax.set_title(
        (f"Ramachandran Plot  |  {title}\n" if title else "Ramachandran Plot\n") +
        f"{n_total} residues  ·  {n_outlier} outliers  ·  {pct_ok:.1f}% favoured+allowed",
        fontsize=10,
    )
    ax.legend(fontsize=8, markerscale=1.5, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved Ramachandran   → {output_path}")
    return output_path


# ─── Buried Surface Area ──────────────────────────────────────────────────────

def _chain_substructure(source: Structure, chain_ids: set[str]) -> Structure:
    """Return a new Structure containing only the specified chain IDs from model 0."""
    from Bio.PDB import Structure as _S, Model as _M
    ns = _S.Structure("sub")
    m  = _M.Model(0)
    ns.add(m)
    for ch in source[0]:
        if ch.id in chain_ids:
            m.add(ch.copy())
    return ns


def buried_surface_areas(
    structure: Structure,
    probe_radius: float = 1.4,
    n_points: int = 100,
) -> list[dict]:
    """
    Compute buried surface area (BSA) for every pair of chains.

    BSA(A, B) = ( SASA(A) + SASA(B) − SASA(A+B) ) / 2

    Also reports how much surface is buried on each individual chain, and
    per-chain interface percentage (BSA_on_X / SASA_X × 100).

    Parameters
    ----------
    structure :
        BioPython Structure (first model used).
    probe_radius :
        Rolling sphere radius in Å. Default 1.4 Å (water).
    n_points :
        Number of sphere points for Shrake-Rupley. Higher = more accurate.
        Default 100.

    Returns
    -------
    List of dicts (one per chain pair), sorted by bsa_total descending:
      chain_1, chain_2,
      sasa_1, sasa_2, sasa_complex,
      bsa_total, bsa_on_1, bsa_on_2,
      interface_pct_1, interface_pct_2
    """
    from itertools import combinations
    from Bio.PDB.SASA import ShrakeRupley

    sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_points)
    chains = [ch.id for ch in structure[0]]

    if len(chains) < 2:
        return []

    # Pre-compute individual chain SASAs
    sasa_single: dict[str, float] = {}
    for cid in chains:
        sub = _chain_substructure(structure, {cid})
        sr.compute(sub, level="A")
        sasa_single[cid] = sum(a.sasa for a in sub.get_atoms())

    results: list[dict] = []
    for ca, cb in combinations(chains, 2):
        sub_ab  = _chain_substructure(structure, {ca, cb})
        sr.compute(sub_ab, level="A")

        sasa_ab = sum(a.sasa for a in sub_ab.get_atoms())
        sasa_a_in_complex = sum(
            a.sasa for ch in sub_ab[0] if ch.id == ca for a in ch.get_atoms()
        )
        sasa_b_in_complex = sasa_ab - sasa_a_in_complex

        sa, sb    = sasa_single[ca], sasa_single[cb]
        bsa_total = (sa + sb - sasa_ab) / 2
        bsa_on_a  = sa - sasa_a_in_complex
        bsa_on_b  = sb - sasa_b_in_complex

        results.append({
            "chain_1":         ca,
            "chain_2":         cb,
            "sasa_1":          round(sa,        1),
            "sasa_2":          round(sb,        1),
            "sasa_complex":    round(sasa_ab,   1),
            "bsa_total":       round(bsa_total, 1),
            "bsa_on_1":        round(bsa_on_a,  1),
            "bsa_on_2":        round(bsa_on_b,  1),
            "interface_pct_1": round(100 * bsa_on_a / sa if sa > 0 else 0.0, 1),
            "interface_pct_2": round(100 * bsa_on_b / sb if sb > 0 else 0.0, 1),
        })

    return sorted(results, key=lambda r: r["bsa_total"], reverse=True)


def plot_bsa_matrix(
    bsa_results: list[dict],
    chain_ids: list[str],
    output_path: str | Path,
    title: str = "",
) -> Path:
    """
    Save a BSA matrix heatmap (chains × chains) and a ranked bar chart.

    Parameters
    ----------
    bsa_results :
        Output of buried_surface_areas().
    chain_ids :
        Ordered list of chain IDs for axis labels.
    output_path :
        PNG output path.
    title :
        Plot title prefix.

    Returns
    -------
    Path to the saved image.
    """
    import numpy as np

    n = len(chain_ids)
    idx = {c: i for i, c in enumerate(chain_ids)}
    matrix = np.zeros((n, n))

    for r in bsa_results:
        i, j = idx[r["chain_1"]], idx[r["chain_2"]]
        matrix[i, j] = r["bsa_total"]
        matrix[j, i] = r["bsa_total"]   # symmetric

    fig_w = min(60, max(8, n * 1.4 + 4))
    fig_h = min(40, max(5, n * 1.2))
    dpi_bsa = min(150, int((150_000_000 / (fig_w * fig_h)) ** 0.5))
    dpi_bsa = max(72, dpi_bsa)
    fig, (ax_heat, ax_bar) = plt.subplots(
        1, 2, figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [n, max(2, n // 2)]},
    )

    # ── Heatmap ──────────────────────────────────────────────────────────────
    vmax = matrix.max() if matrix.max() > 0 else 1.0
    im = ax_heat.imshow(matrix, cmap="YlOrRd", vmin=0, vmax=vmax, aspect="equal")
    plt.colorbar(im, ax=ax_heat, label="BSA (Å²)", shrink=0.8)

    ax_heat.set_xticks(range(n))
    ax_heat.set_yticks(range(n))
    ax_heat.set_xticklabels([f"Chain {c}" for c in chain_ids], fontsize=9)
    ax_heat.set_yticklabels([f"Chain {c}" for c in chain_ids], fontsize=9)

    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            if val > 0:
                txt_color = "white" if val > 0.6 * vmax else "black"
                ax_heat.text(j, i, f"{val:.0f}", ha="center", va="center",
                             fontsize=8, color=txt_color, fontweight="bold")

    ax_heat.set_title(
        (f"{title}  —  " if title else "") + "Buried Surface Area (Å²)",
        fontsize=11, fontweight="bold",
    )

    # ── Ranked bar chart ──────────────────────────────────────────────────────
    ranked = [r for r in bsa_results if r["bsa_total"] > 0]
    if ranked:
        labels = [f"{r['chain_1']}–{r['chain_2']}" for r in ranked]
        values = [r["bsa_total"]                   for r in ranked]
        colors = plt.cm.YlOrRd(
            [v / max(values) * 0.85 + 0.1 for v in values]
        )
        bars = ax_bar.barh(labels[::-1], values[::-1],
                           color=colors[::-1], edgecolor="white", linewidth=0.6)
        ax_bar.bar_label(bars, fmt="%.0f", fontsize=8, padding=3)
        ax_bar.set_xlabel("BSA (Å²)", fontsize=9)
        ax_bar.set_title("Interface ranking", fontsize=10, fontweight="bold")
        ax_bar.set_xlim(0, max(values) * 1.25)
        ax_bar.spines[["top", "right"]].set_visible(False)
        ax_bar.tick_params(axis="y", labelsize=9)
    else:
        ax_bar.text(0.5, 0.5, "No buried interfaces\ndetected",
                    ha="center", va="center", transform=ax_bar.transAxes,
                    fontsize=11, color="#888")
        ax_bar.set_axis_off()

    plt.tight_layout()
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=dpi_bsa, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved BSA matrix   → {output_path}")
    return output_path


# ─── Batch Summary ────────────────────────────────────────────────────────────

def batch_summary(paths: list[str | Path],
                   out_csv: str | Path | None = None) -> list[dict]:
    """
    Summarize a list of structure files and optionally save as CSV.

    Parameters
    ----------
    paths :
        List of PDB or CIF file paths.
    out_csv :
        If provided, write results to this CSV file.

    Returns
    -------
    List of summary dicts (one per structure).
    """
    rows: list[dict] = []
    for path in paths:
        path = Path(path)
        print(f"  Summarizing {path.name} …", end=" ", flush=True)
        try:
            row = summarize_structure(path)
            rows.append(row)
            print("OK")
        except Exception as e:
            print(f"ERROR: {e}")

    if out_csv and rows:
        out_csv = Path(out_csv)
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"\n  Saved batch summary → {out_csv}  ({len(rows)} entries)")

    return rows


def print_summary(row: dict) -> None:
    """Print a single summary dict to the terminal."""
    print(f"\n{'─'*55}")
    print(f"  PDB:        {row['pdb_id']}")
    print(f"  Title:      {row['title'][:60]}")
    print(f"  Method:     {row['method']}")
    print(f"  Resolution: {row['resolution']} Å")
    print(f"  R-work:     {row['r_work']}   R-free: {row['r_free']}")
    print(f"  Organism:   {row['organism']}")
    print(f"  Deposited:  {row['deposition_date']}")
    print(f"  Models:     {row['n_models']}")
    print(f"  Chains:     {row['chains']}  ({row['residues_per_chain']})")
    print(f"  Residues:   {row['n_residues']}   Atoms: {row['n_atoms']}")
    print(f"  Ligands:    {row['n_ligands']}   ({row['ligand_names']})")
    print(f"  Waters:     {row['n_waters']}")
    miss = row['n_missing_residues']
    print(f"  Missing:    {miss} chain(s) with gaps" if miss else "  Missing:    none detected")
    print(f"{'─'*55}")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Structure inventory, metadata, geometry, and batch summaries.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 summarize_structures.py 1abc.pdb
  python3 summarize_structures.py Files/*.pdb --out summary.csv
  python3 summarize_structures.py 1abc.pdb --missing
  python3 summarize_structures.py 1abc.pdb --fasta
  python3 summarize_structures.py 1abc.pdb --fasta --chain A
  python3 summarize_structures.py 1abc.pdb --bfactor --plot
  python3 summarize_structures.py 1abc.pdb --ramachandran --plot
  python3 summarize_structures.py 1abc.pdb --rg
""",
    )

    parser.add_argument("files", nargs="+", help="PDB or CIF files")

    modes = parser.add_argument_group("mode (default: print summary)")
    modes.add_argument("--missing",      action="store_true",
                       help="Show missing residue report")
    modes.add_argument("--fasta",        action="store_true",
                       help="Print FASTA sequence from ATOM records")
    modes.add_argument("--bfactor",      action="store_true",
                       help="Print per-residue B-factor table")
    modes.add_argument("--ramachandran", action="store_true",
                       help="Print phi/psi table with outlier flags")
    modes.add_argument("--rg",           action="store_true",
                       help="Print radius of gyration")

    opts = parser.add_argument_group("options")
    opts.add_argument("--chain", default=None,
                      help="Restrict analysis to one chain")
    opts.add_argument("--plot",  action="store_true",
                      help="Save PNG plot for --bfactor or --ramachandran")
    opts.add_argument("--out",   default=None, metavar="CSV",
                      help="Save batch summary to CSV (default mode only)")
    opts.add_argument("--outdir", default=".", metavar="DIR",
                      help="Directory for plot output (default: .)")

    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    paths = [Path(p) for p in args.files]

    # ── FASTA ─────────────────────────────────────────────────────────────────
    if args.fasta:
        for path in paths:
            print(extract_fasta(path, args.chain))
        return

    # ── Missing residues ──────────────────────────────────────────────────────
    if args.missing:
        for path in paths:
            print(f"\n{path.name}")
            rows = find_missing_residues(path)
            if not rows:
                print("  No SEQRES records found or no gaps detected.")
            else:
                print(f"  {'Chain':<6} {'Expected':>10} {'Observed':>10} "
                      f"{'Missing':>9} {'%Modeled':>10}")
                print("  " + "─" * 50)
                for r in rows:
                    print(f"  {r['chain']:<6} {r['expected_residues']:>10} "
                          f"{r['observed_residues']:>10} {r['missing_count']:>9} "
                          f"{r['pct_modeled']:>9.1f}%")
        return

    # ── B-factor ──────────────────────────────────────────────────────────────
    if args.bfactor:
        for path in paths:
            structure = load_structure(path)
            data = bfactor_stats(structure, args.chain)
            print(f"\n{path.name}  ({len(data)} residues)")
            print(f"  {'Chain':<6} {'Resi':>5} {'Resn':<5} {'B-factor':>9} {'Flag'}")
            print("  " + "─" * 36)
            for r in data:
                flag = f"  ← {r['flag']}" if r["flag"] else ""
                print(f"  {r['chain']:<6} {r['resi']:>5} {r['resn']:<5} "
                      f"{r['bfactor']:>9.2f}{flag}")
            if args.plot:
                stem = path.stem
                plot_bfactor(data, outdir / f"{stem}_bfactor.png",
                             title=path.stem.upper())
        return

    # ── Ramachandran ──────────────────────────────────────────────────────────
    if args.ramachandran:
        for path in paths:
            structure = load_structure(path)
            data = ramachandran_data(structure, args.chain)
            outliers = [r for r in data if r["region"] == "outlier"]
            pct = 100 * (1 - len(outliers) / len(data)) if data else 0
            print(f"\n{path.name}  {len(data)} residues  "
                  f"{len(outliers)} outliers  {pct:.1f}% favored+allowed")
            if outliers:
                print(f"\n  Outliers:")
                print(f"  {'Chain':<6} {'Resi':>5} {'Resn':<5} {'phi':>8} {'psi':>8}")
                print("  " + "─" * 36)
                for r in outliers:
                    print(f"  {r['chain']:<6} {r['resi']:>5} {r['resn']:<5} "
                          f"{r['phi']:>8.1f} {r['psi']:>8.1f}")
            if args.plot:
                plot_ramachandran(data, outdir / f"{path.stem}_rama.png",
                                  title=path.stem.upper())
        return

    # ── Radius of gyration ────────────────────────────────────────────────────
    if args.rg:
        for path in paths:
            structure = load_structure(path)
            rg = radius_of_gyration(structure, args.chain)
            print(f"{path.name}  Rg = {rg:.2f} Å")
        return

    # ── Default: summary (single or batch) ───────────────────────────────────
    if len(paths) == 1:
        row = summarize_structure(paths[0])
        print_summary(row)
        if args.out:
            batch_summary(paths, args.out)
    else:
        batch_summary(paths, args.out)


if __name__ == "__main__":
    main()
