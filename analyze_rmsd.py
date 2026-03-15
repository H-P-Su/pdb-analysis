#!/usr/bin/env python3
"""
analyze_rmsd.py — RMSD calculation, structure superposition, and ensemble analysis.

All functions importable for use in other projects.

Functions
---------
load_ca_atoms(path, chain_id)           Load Cα atoms + sequence string
align_sequences(seq1, seq2)             Pairwise sequence alignment → index pairs
calculate_rmsd(path1, path2, ...)       Cα RMSD between two structures (auto-aligns)
superpose(path1, path2, ...)            Superpose struct2 onto struct1, save optional
per_residue_rmsd(path1, path2, ...)     Per-residue Cα distances after superposition
pairwise_rmsd_matrix(paths, ...)        All-vs-all RMSD matrix + labels
nmr_ensemble_rmsd(path, ...)            Inter-model RMSD for NMR multi-MODEL files
plot_rmsd_matrix(matrix, labels, ...)   Heatmap → PNG
plot_per_residue_rmsd(values, ...)      Line plot → PNG

CLI
---
python3 analyze_rmsd.py struct1.pdb struct2.pdb
python3 analyze_rmsd.py struct1.pdb struct2.pdb --chain A --per-residue --plot
python3 analyze_rmsd.py struct1.pdb struct2.pdb --superpose --out superposed.pdb
python3 analyze_rmsd.py Files/*.pdb --matrix --out rmsd_matrix.csv --plot
python3 analyze_rmsd.py nmr_ensemble.pdb --ensemble
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Superimposer
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure

# ─── Constants ────────────────────────────────────────────────────────────────

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "SEC": "U",
}

PROTEIN_RESIDUES = frozenset(THREE_TO_ONE.keys())


# ─── Structure Loading ─────────────────────────────────────────────────────────

def _load(path: str | Path, sid: str = "s") -> Structure:
    path = Path(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True).get_structure(sid, str(path))
    return PDBParser(QUIET=True).get_structure(sid, str(path))


# ─── Cα Extraction ────────────────────────────────────────────────────────────

def load_ca_atoms(path: str | Path,
                   chain_id: str | None = None,
                   model_id: int = 0) -> tuple[list[Atom], str]:
    """
    Load Cα atoms and one-letter sequence from a structure file.

    Parameters
    ----------
    path :
        PDB or CIF file.
    chain_id :
        Restrict to one chain.  If None, concatenates all chains.
    model_id :
        Model index (0-based).  Useful for NMR ensembles.

    Returns
    -------
    (atoms, sequence) where atoms[i] corresponds to sequence[i].
    """
    structure = _load(path, Path(path).stem)
    models = list(structure.get_models())
    if model_id >= len(models):
        raise ValueError(f"Model {model_id} not found (structure has {len(models)} models)")
    model = models[model_id]

    atoms: list[Atom] = []
    seq_chars: list[str] = []

    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for res in chain:
            if res.id[0] != " ":
                continue
            resn = res.get_resname().strip()
            if resn not in PROTEIN_RESIDUES:
                continue
            if "CA" not in res:
                continue
            atoms.append(res["CA"])
            seq_chars.append(THREE_TO_ONE.get(resn, "X"))

    return atoms, "".join(seq_chars)


# ─── Sequence Alignment ────────────────────────────────────────────────────────

def align_sequences(seq1: str, seq2: str) -> tuple[list[int], list[int]]:
    """
    Global pairwise sequence alignment; return matching index pairs.

    Parameters
    ----------
    seq1, seq2 :
        One-letter amino-acid sequences.

    Returns
    -------
    (idx1, idx2) — parallel lists of matched positions in seq1 and seq2.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = 0
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = list(aligner.align(seq1, seq2))
    if not alignments:
        return [], []

    best = alignments[0]
    idx1: list[int] = []
    idx2: list[int] = []

    # aligned attribute: tuple of two tuples of (start, end) block pairs
    for (s1, e1), (s2, e2) in zip(*best.aligned):
        for i, j in zip(range(s1, e1), range(s2, e2)):
            idx1.append(i)
            idx2.append(j)

    return idx1, idx2


# ─── RMSD ─────────────────────────────────────────────────────────────────────

def calculate_rmsd(path1: str | Path,
                    path2: str | Path,
                    chain1: str | None = None,
                    chain2: str | None = None,
                    align_by: str = "sequence",
                    model1: int = 0,
                    model2: int = 0) -> float:
    """
    Compute Cα RMSD between two structures after optimal superposition.

    Parameters
    ----------
    path1, path2 :
        PDB or CIF files.
    chain1, chain2 :
        Chain IDs to use (None = all chains concatenated).
    align_by :
        ``"sequence"`` — global pairwise alignment (default, handles
        insertions/deletions).
        ``"resnum"``   — match by residue number (faster; requires same numbering).
    model1, model2 :
        Model indices (for NMR ensembles).

    Returns
    -------
    Cα RMSD in Å.
    """
    atoms1, seq1 = load_ca_atoms(path1, chain1, model1)
    atoms2, seq2 = load_ca_atoms(path2, chain2, model2)

    fixed, moving = _match_atoms(atoms1, seq1, atoms2, seq2, align_by, path1, path2)

    if len(fixed) < 3:
        raise ValueError(
            f"Only {len(fixed)} matched Cα atoms — too few to compute RMSD.")

    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    return round(float(sup.rms), 4)


def _match_atoms(atoms1: list[Atom], seq1: str,
                  atoms2: list[Atom], seq2: str,
                  align_by: str,
                  path1: str | Path,
                  path2: str | Path) -> tuple[list[Atom], list[Atom]]:
    """Return parallel lists of matched Cα atoms."""
    if align_by == "resnum":
        # Match by residue sequence number
        rn1 = {a.get_parent().id[1]: a for a in atoms1}
        rn2 = {a.get_parent().id[1]: a for a in atoms2}
        common = sorted(set(rn1) & set(rn2))
        return [rn1[r] for r in common], [rn2[r] for r in common]

    # sequence alignment
    idx1, idx2 = align_sequences(seq1, seq2)
    if not idx1:
        # fallback: pair by position
        n = min(len(atoms1), len(atoms2))
        return atoms1[:n], atoms2[:n]
    return [atoms1[i] for i in idx1], [atoms2[i] for i in idx2]


# ─── Superposition ────────────────────────────────────────────────────────────

def superpose(path1: str | Path,
               path2: str | Path,
               chain1: str | None = None,
               chain2: str | None = None,
               align_by: str = "sequence",
               output_path: str | Path | None = None) -> tuple[float, Structure]:
    """
    Superpose structure 2 onto structure 1 (minimise Cα RMSD).

    Parameters
    ----------
    path1, path2 :
        Reference and mobile PDB/CIF files.
    chain1, chain2 :
        Chain IDs (None = all).
    align_by :
        ``"sequence"`` or ``"resnum"``.
    output_path :
        If given, save the superposed structure 2 to this PDB file.

    Returns
    -------
    (rmsd, superposed_structure) — RMSD in Å and the transformed BioPython Structure.
    """
    atoms1, seq1 = load_ca_atoms(path1, chain1)
    atoms2, seq2 = load_ca_atoms(path2, chain2)

    fixed, moving = _match_atoms(atoms1, seq1, atoms2, seq2, align_by, path1, path2)

    sup = Superimposer()
    sup.set_atoms(fixed, moving)

    structure2 = _load(path2, "mobile")
    sup.apply(list(structure2.get_atoms()))

    rmsd = round(float(sup.rms), 4)

    if output_path:
        output_path = Path(output_path)
        io = PDBIO()
        io.set_structure(structure2)
        io.save(str(output_path))
        print(f"  Superposed structure → {output_path}  (RMSD {rmsd} Å)")

    return rmsd, structure2


# ─── Per-residue RMSD ────────────────────────────────────────────────────────

def per_residue_rmsd(path1: str | Path,
                      path2: str | Path,
                      chain1: str | None = None,
                      chain2: str | None = None,
                      align_by: str = "sequence") -> list[dict]:
    """
    Compute per-residue Cα distance after global superposition.

    Parameters
    ----------
    path1, path2 :
        PDB/CIF files.
    chain1, chain2 :
        Chain IDs (None = all).
    align_by :
        ``"sequence"`` or ``"resnum"``.

    Returns
    -------
    List of dicts: {resi, resn, chain, rmsd_ca}
    """
    atoms1, seq1 = load_ca_atoms(path1, chain1)
    atoms2, seq2 = load_ca_atoms(path2, chain2)

    fixed, moving = _match_atoms(atoms1, seq1, atoms2, seq2, align_by, path1, path2)

    sup = Superimposer()
    sup.set_atoms(fixed, moving)

    # Apply transformation to moving atoms
    # Superimposer.apply works on a list; we only want to move the matched atoms
    # Create a temporary copy via numpy
    rot, tran = sup.rotran
    rows: list[dict] = []
    for a_fix, a_mov in zip(fixed, moving):
        coords_mov = a_mov.get_vector().get_array()
        coords_sup = np.dot(rot, coords_mov) + tran
        coords_fix = a_fix.get_vector().get_array()
        dist = float(np.linalg.norm(coords_fix - coords_sup))
        res = a_fix.get_parent()
        rows.append({
            "chain":   res.get_parent().id,
            "resi":    res.id[1],
            "resn":    res.get_resname().strip(),
            "rmsd_ca": round(dist, 3),
        })

    return rows


def plot_per_residue_rmsd(data: list[dict],
                           output_path: str | Path,
                           title: str = "",
                           threshold: float = 2.0) -> Path:
    """
    Save a per-residue RMSD line plot.

    Parameters
    ----------
    data :
        Output of per_residue_rmsd().
    output_path :
        PNG output path.
    title :
        Plot title.
    threshold :
        Mark residues above this RMSD (Å) in red. Default 2.0 Å.
    """
    if not data:
        print("  No per-residue RMSD data to plot.")
        return Path(output_path)

    x     = [r["resi"]    for r in data]
    y     = [r["rmsd_ca"] for r in data]
    resns = [r["resn"]    for r in data]

    colors = ["#e84040" if v > threshold else "#4488cc" for v in y]

    fig, ax = plt.subplots(figsize=(max(10, len(x) * 0.07 + 2), 4))
    ax.bar(x, y, color=colors, width=1.0, linewidth=0)
    ax.axhline(threshold, color="#e84040", linewidth=1.0, linestyle="--",
               label=f"Threshold {threshold} Å")
    ax.axhline(float(np.mean(y)), color="black", linewidth=0.8, linestyle=":",
               label=f"Mean {np.mean(y):.2f} Å")

    # Label high-RMSD residues
    high = [(xi, yi, rn) for xi, yi, rn in zip(x, y, resns) if yi > threshold]
    for xi, yi, rn in high[:15]:
        ax.annotate(f"{rn}{xi}", (xi, yi),
                    textcoords="offset points", xytext=(0, 4),
                    fontsize=6, ha="center", color="#e84040")

    ax.set_xlabel("Residue number", fontsize=10)
    ax.set_ylabel("Cα RMSD (Å)", fontsize=10)
    ax.set_title((f"{title}\n" if title else "") +
                 f"Per-residue Cα RMSD after superposition  "
                 f"({sum(1 for v in y if v > threshold)} residues > {threshold} Å)",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved per-residue RMSD → {output_path}")
    return output_path


# ─── Pairwise RMSD Matrix ─────────────────────────────────────────────────────

def pairwise_rmsd_matrix(paths: list[str | Path],
                          chain: str | None = None,
                          align_by: str = "sequence") -> tuple[np.ndarray, list[str]]:
    """
    Compute all-vs-all Cα RMSD matrix for a set of structures.

    Parameters
    ----------
    paths :
        List of PDB/CIF files.
    chain :
        Chain ID to use (None = all chains).
    align_by :
        ``"sequence"`` or ``"resnum"``.

    Returns
    -------
    (matrix, labels) — symmetric n×n RMSD matrix and structure labels.
    """
    n = len(paths)
    labels = [Path(p).stem for p in paths]
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            try:
                rmsd = calculate_rmsd(paths[i], paths[j],
                                       chain1=chain, chain2=chain,
                                       align_by=align_by)
            except Exception as e:
                print(f"  RMSD({labels[i]}, {labels[j]}): {e}")
                rmsd = float("nan")
            matrix[i, j] = rmsd
            matrix[j, i] = rmsd
            print(f"  RMSD  {labels[i]:20s} vs {labels[j]:20s}  →  {rmsd:.3f} Å")

    return matrix, labels


def plot_rmsd_matrix(matrix: np.ndarray,
                      labels: list[str],
                      output_path: str | Path,
                      title: str = "Pairwise Cα RMSD (Å)") -> Path:
    """
    Save a heatmap of the pairwise RMSD matrix.

    Parameters
    ----------
    matrix :
        n×n symmetric RMSD matrix.
    labels :
        Structure labels (length n).
    output_path :
        PNG output path.
    title :
        Plot title.
    """
    n = len(labels)
    fig_size = max(5, n * 0.7 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    # Mask diagonal for cleaner colour scaling
    display = matrix.copy()
    np.fill_diagonal(display, np.nan)

    im = ax.imshow(display, cmap="coolwarm", interpolation="nearest",
                   vmin=0, vmax=np.nanmax(display) or 1)
    plt.colorbar(im, ax=ax, shrink=0.8, label="RMSD (Å)")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(labels, fontsize=9)

    for i in range(n):
        for j in range(n):
            if i != j and not np.isnan(matrix[i, j]):
                ax.text(j, i, f"{matrix[i, j]:.2f}",
                        ha="center", va="center", fontsize=8,
                        color="white" if matrix[i, j] > np.nanmax(display) * 0.6 else "black")

    ax.set_title(title, fontsize=12, pad=10)
    plt.tight_layout()
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved RMSD matrix    → {output_path}")
    return output_path


# ─── NMR Ensemble ─────────────────────────────────────────────────────────────

def nmr_ensemble_rmsd(path: str | Path,
                       chain_id: str | None = None,
                       reference_model: int = 0) -> list[dict]:
    """
    Compute per-model Cα RMSD in an NMR multi-MODEL ensemble.

    Each model is superposed onto the reference model (default: model 0).

    Parameters
    ----------
    path :
        PDB file with multiple MODEL records.
    chain_id :
        Chain to analyse (None = all).
    reference_model :
        Index of the reference model (0-based).

    Returns
    -------
    List of dicts: {model_id, n_atoms, rmsd}
    """
    structure = _load(path, "nmr")
    models = list(structure.get_models())
    n = len(models)
    if n < 2:
        print("  Only one model found — not an NMR ensemble.")
        return []

    ref_atoms, ref_seq = load_ca_atoms(path, chain_id, reference_model)
    rows: list[dict] = []

    for idx in range(n):
        if idx == reference_model:
            rows.append({"model_id": idx, "n_atoms": len(ref_atoms), "rmsd": 0.0})
            continue
        try:
            mob_atoms, mob_seq = load_ca_atoms(path, chain_id, idx)
            fixed, moving = _match_atoms(ref_atoms, ref_seq,
                                          mob_atoms, mob_seq, "sequence", path, path)
            sup = Superimposer()
            sup.set_atoms(fixed, moving)
            rows.append({"model_id": idx, "n_atoms": len(fixed),
                         "rmsd": round(float(sup.rms), 4)})
        except Exception as e:
            rows.append({"model_id": idx, "n_atoms": 0, "rmsd": float("nan")})
            print(f"  Model {idx}: {e}")

    return rows


def plot_ensemble_rmsd(data: list[dict],
                        output_path: str | Path,
                        title: str = "") -> Path:
    """
    Bar chart of per-model RMSD for an NMR ensemble.

    Parameters
    ----------
    data :
        Output of nmr_ensemble_rmsd().
    output_path :
        PNG output path.
    title :
        Plot title (optional).
    """
    if not data:
        return Path(output_path)

    models = [r["model_id"] for r in data]
    rmsds  = [r["rmsd"] for r in data]
    mean_r = float(np.nanmean([r for r in rmsds if r > 0]))

    fig, ax = plt.subplots(figsize=(max(6, len(models) * 0.4 + 2), 4))
    colors = ["#cccccc" if r == 0 else
              "#e84040" if r > mean_r + np.nanstd(rmsds) else "#4488cc"
              for r in rmsds]
    ax.bar(models, rmsds, color=colors, edgecolor="white", linewidth=0.5)
    ax.axhline(mean_r, color="black", linewidth=1.0, linestyle="--",
               label=f"Mean {mean_r:.2f} Å")

    ax.set_xlabel("Model index", fontsize=10)
    ax.set_ylabel("Cα RMSD vs. model 0 (Å)", fontsize=10)
    ax.set_title((f"{title}\n" if title else "") +
                 f"NMR Ensemble RMSD  ({len(models)} models)",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved ensemble RMSD  → {output_path}")
    return output_path


# ─── CSV Export ───────────────────────────────────────────────────────────────

def save_csv(rows: list[dict], path: str | Path) -> None:
    path = Path(path)
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved CSV            → {path}")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="RMSD calculation, superposition, and ensemble analysis.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # RMSD between two structures
  python3 analyze_rmsd.py struct1.pdb struct2.pdb

  # RMSD restricted to chain A, with per-residue plot
  python3 analyze_rmsd.py struct1.pdb struct2.pdb --chain A --per-residue --plot

  # Superpose struct2 onto struct1, save result
  python3 analyze_rmsd.py struct1.pdb struct2.pdb --superpose --out superposed.pdb

  # All-vs-all matrix for a set of structures
  python3 analyze_rmsd.py Files/*.pdb --matrix --out rmsd.csv --plot

  # NMR ensemble RMSD
  python3 analyze_rmsd.py nmr.pdb --ensemble --plot
""",
    )

    parser.add_argument("files", nargs="+", help="PDB or CIF files")

    modes = parser.add_argument_group("mode")
    mx = modes.add_mutually_exclusive_group()
    mx.add_argument("--matrix",      action="store_true",
                    help="All-vs-all pairwise RMSD matrix (≥2 files)")
    mx.add_argument("--per-residue", action="store_true", dest="per_residue",
                    help="Per-residue RMSD after superposition (exactly 2 files)")
    mx.add_argument("--superpose",   action="store_true",
                    help="Superpose file2 onto file1 (exactly 2 files)")
    mx.add_argument("--ensemble",    action="store_true",
                    help="NMR ensemble RMSD (1 multi-MODEL file)")

    opts = parser.add_argument_group("options")
    opts.add_argument("--chain",  default=None, help="Chain ID to use")
    opts.add_argument("--chain2", default=None,
                      help="Chain ID for second structure (default: same as --chain)")
    opts.add_argument("--align",  default="sequence",
                      choices=["sequence", "resnum"],
                      help="Alignment method (default: sequence)")
    opts.add_argument("--threshold", type=float, default=2.0,
                      help="RMSD threshold for per-residue highlighting (Å, default 2.0)")
    opts.add_argument("--ref-model", type=int, default=0, dest="ref_model",
                      help="Reference model index for NMR ensemble (default: 0)")

    out = parser.add_argument_group("output")
    out.add_argument("--out",    default=None,
                     help="CSV output (matrix/per-residue) or PDB (superpose)")
    out.add_argument("--plot",   action="store_true",
                     help="Save PNG plot")
    out.add_argument("--outdir", default=".", metavar="DIR",
                     help="Directory for plots (default: .)")

    args = parser.parse_args()
    paths  = [Path(p) for p in args.files]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chain2 = args.chain2 or args.chain

    # ── NMR ensemble ──────────────────────────────────────────────────────────
    if args.ensemble:
        if len(paths) != 1:
            parser.error("--ensemble requires exactly one file")
        data = nmr_ensemble_rmsd(paths[0], args.chain, args.ref_model)
        if not data:
            return
        rmsds = [r["rmsd"] for r in data if r["rmsd"] > 0]
        print(f"\n{'Model':>6}  {'Atoms':>6}  {'RMSD (Å)':>10}")
        print("─" * 30)
        for r in data:
            print(f"{r['model_id']:>6}  {r['n_atoms']:>6}  {r['rmsd']:>10.4f}")
        if rmsds:
            print(f"\n  Mean RMSD: {np.mean(rmsds):.3f} Å  "
                  f"Max: {np.max(rmsds):.3f} Å")
        if args.out:
            save_csv(data, args.out)
        if args.plot:
            plot_ensemble_rmsd(data, outdir / f"{paths[0].stem}_ensemble_rmsd.png",
                               title=paths[0].stem.upper())
        return

    # ── Pairwise matrix ────────────────────────────────────────────────────────
    if args.matrix:
        if len(paths) < 2:
            parser.error("--matrix requires at least 2 files")
        matrix, labels = pairwise_rmsd_matrix(paths, args.chain, args.align)
        print(f"\nPairwise Cα RMSD matrix (Å):\n")
        w = max(len(l) for l in labels) + 2
        print(" " * w + "  ".join(f"{l:>{w}}" for l in labels))
        for i, row_label in enumerate(labels):
            vals = "  ".join(
                "  -   " if j == i else f"{matrix[i,j]:>{w}.3f}"
                for j in range(len(labels))
            )
            print(f"{row_label:>{w}}  {vals}")

        if args.out:
            csv_path = Path(args.out)
            with open(csv_path, "w", newline="") as f:
                w_csv = csv.writer(f)
                w_csv.writerow([""] + labels)
                for i, label in enumerate(labels):
                    w_csv.writerow([label] + [f"{matrix[i,j]:.4f}"
                                               for j in range(len(labels))])
            print(f"\n  Saved matrix CSV     → {csv_path}")
        if args.plot:
            plot_rmsd_matrix(matrix, labels,
                              outdir / "rmsd_matrix.png")
        return

    # ── Two-structure modes ────────────────────────────────────────────────────
    if len(paths) < 2:
        parser.error("Provide at least 2 files (or use --ensemble for NMR)")

    path1, path2 = paths[0], paths[1]

    if args.superpose:
        out_pdb = Path(args.out) if args.out else outdir / f"{path2.stem}_superposed.pdb"
        rmsd, _ = superpose(path1, path2, args.chain, chain2,
                             args.align, out_pdb)
        print(f"\n  Cα RMSD after superposition: {rmsd:.4f} Å")
        return

    if args.per_residue:
        data = per_residue_rmsd(path1, path2, args.chain, chain2, args.align)
        print(f"\n{'Chain':<6} {'Resi':>5} {'Resn':<5} {'RMSD_Cα (Å)':>12}")
        print("─" * 32)
        for r in data:
            flag = "  ←" if r["rmsd_ca"] > args.threshold else ""
            print(f"{r['chain']:<6} {r['resi']:>5} {r['resn']:<5} "
                  f"{r['rmsd_ca']:>12.3f}{flag}")
        overall = float(np.sqrt(np.mean([r["rmsd_ca"]**2 for r in data])))
        print(f"\n  Overall RMSD: {overall:.4f} Å  ({len(data)} aligned residues)")
        if args.out:
            save_csv(data, args.out)
        if args.plot:
            title = f"{path1.stem} vs {path2.stem}"
            plot_per_residue_rmsd(data, outdir / f"{path1.stem}_vs_{path2.stem}_rmsd.png",
                                  title=title, threshold=args.threshold)
        return

    # ── Default: single RMSD value ─────────────────────────────────────────────
    rmsd = calculate_rmsd(path1, path2, args.chain, chain2, args.align)
    n1 = len(load_ca_atoms(path1, args.chain)[0])
    n2 = len(load_ca_atoms(path2, chain2)[0])
    print(f"\n  {path1.name}  ({n1} Cα)")
    print(f"  {path2.name}  ({n2} Cα)")
    print(f"\n  Cα RMSD: {rmsd:.4f} Å  (alignment: {args.align})")


if __name__ == "__main__":
    main()
