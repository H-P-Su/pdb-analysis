#!/usr/bin/env python3
"""
ligand_rmsd.py — Compute ligand RMSD across multiple PDB structures.

Aligns each structure to a reference via protein backbone, then computes
heavy-atom RMSD for a specified ligand residue. Outliers are flagged by
z-score and excluded from the final average (threshold is configurable).

If no reference is provided, a centroid structure (mean coordinates across
all structures) is used.

Usage:
    python ligand_rmsd.py --ligand LIG "structures/*.pdb"
    python ligand_rmsd.py --ligand LIG --reference ref.pdb "structures/*.pdb"
    python ligand_rmsd.py --ligand LIG --chain A --zscore 1.5 "structures/*.pdb"

Dependencies:
    pip install biopython numpy tabulate
"""

import argparse
import glob
import sys
import warnings
from pathlib import Path

import numpy as np
from tabulate import tabulate

# Suppress Bio warnings about discontinuous chains etc.
warnings.filterwarnings("ignore")

from Bio import PDB
from Bio.PDB import Superimposer


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def load_structure(path: str, parser: PDB.PDBParser) -> PDB.Structure.Structure:
    struct_id = Path(path).stem
    return parser.get_structure(struct_id, path)


def get_backbone_atoms(structure, chain_id: str | None = None) -> list:
    """Return CA atoms from the first model, optionally filtered by chain."""
    backbone = []
    model = structure[0]
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for residue in chain:
            if residue.id[0] != " ":  # skip HETATM / water
                continue
            if "CA" in residue:
                backbone.append(residue["CA"])
    return backbone


def get_ligand_heavy_atoms(structure, resn: str, chain_id: str | None = None) -> dict:
    """
    Return dict of {atom_name: np.array([x,y,z])} for heavy atoms of the
    specified ligand residue name. Raises if ligand not found.
    """
    model = structure[0]
    atoms = {}
    found_residue = False

    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for residue in chain:
            if residue.resname.strip() != resn.strip():
                continue
            found_residue = True
            for atom in residue.get_atoms():
                element = atom.element.strip() if atom.element else atom.name[0]
                if element == "H":
                    continue
                atoms[atom.name] = atom.get_vector().get_array()

    if not found_residue:
        raise ValueError(
            f"Ligand '{resn}' not found in structure '{structure.id}'"
            + (f" chain '{chain_id}'" if chain_id else "")
        )
    return atoms


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

def align_to_reference(mobile_struct, ref_struct, chain_id: str | None = None) -> None:
    """Superimpose mobile onto ref in-place using backbone CA atoms."""
    ref_atoms = get_backbone_atoms(ref_struct, chain_id)
    mob_atoms = get_backbone_atoms(mobile_struct, chain_id)

    # Match by residue sequence number + chain
    ref_map = {
        (a.get_parent().get_parent().id, a.get_parent().id[1]): a
        for a in ref_atoms
    }
    mob_map = {
        (a.get_parent().get_parent().id, a.get_parent().id[1]): a
        for a in mob_atoms
    }

    common_keys = sorted(set(ref_map) & set(mob_map))
    if len(common_keys) < 3:
        raise ValueError(
            f"Too few common backbone residues to align "
            f"'{mobile_struct.id}' (found {len(common_keys)})"
        )

    ref_aligned = [ref_map[k] for k in common_keys]
    mob_aligned = [mob_map[k] for k in common_keys]

    sup = Superimposer()
    sup.set_atoms(ref_aligned, mob_aligned)
    sup.apply(mobile_struct.get_atoms())


# ---------------------------------------------------------------------------
# Centroid reference
# ---------------------------------------------------------------------------

def compute_centroid_coords(all_atom_dicts: list[dict]) -> dict:
    """
    Given a list of {atom_name: coords} dicts, return a centroid dict
    containing the mean coords for atoms common to ALL structures.
    """
    common_names = set(all_atom_dicts[0].keys())
    for d in all_atom_dicts[1:]:
        common_names &= set(d.keys())

    centroid = {}
    for name in common_names:
        stack = np.array([d[name] for d in all_atom_dicts])
        centroid[name] = stack.mean(axis=0)
    return centroid


# ---------------------------------------------------------------------------
# RMSD
# ---------------------------------------------------------------------------

def compute_rmsd(coords1: dict, coords2: dict) -> tuple[float, int]:
    """
    Compute RMSD between two atom-name-keyed coord dicts.
    Returns (rmsd, n_atoms_matched).
    """
    common = sorted(set(coords1) & set(coords2))
    if not common:
        raise ValueError("No common atom names found between ligand instances.")

    pts1 = np.array([coords1[n] for n in common])
    pts2 = np.array([coords2[n] for n in common])
    rmsd = float(np.sqrt(np.mean(np.sum((pts1 - pts2) ** 2, axis=1))))
    return rmsd, len(common)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute ligand RMSD across PDB structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "structures",
        nargs="+",
        help="Glob pattern(s) or explicit paths to PDB files.",
    )
    p.add_argument(
        "--ligand", "-l",
        required=True,
        help="Ligand residue name (e.g. LIG, ATP, HEM).",
    )
    p.add_argument(
        "--reference", "-r",
        default=None,
        help="Path to reference PDB. If omitted, centroid of all structures is used.",
    )
    p.add_argument(
        "--chain", "-c",
        default=None,
        help="Protein chain ID to use for alignment (default: all chains).",
    )
    p.add_argument(
        "--zscore", "-z",
        type=float,
        default=2.0,
        help="Z-score threshold for outlier exclusion (default: 2.0).",
    )
    p.add_argument(
        "--no-exclude",
        action="store_true",
        help="Report z-scores but do not exclude any structures.",
    )
    return p.parse_args()


def expand_globs(patterns: list[str]) -> list[str]:
    paths = []
    for pattern in patterns:
        expanded = sorted(glob.glob(pattern))
        if not expanded:
            print(f"Warning: no files matched pattern '{pattern}'", file=sys.stderr)
        paths.extend(expanded)
    # Deduplicate while preserving order
    seen = set()
    result = []
    for p in paths:
        if p not in seen:
            seen.add(p)
            result.append(p)
    return result


def main():
    args = parse_args()
    parser = PDB.PDBParser(QUIET=True)

    # Expand globs
    pdb_files = expand_globs(args.structures)
    if not pdb_files:
        print("Error: no PDB files found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(pdb_files)} structure(s).")

    # Load all structures
    structures = {}
    load_errors = []
    for path in pdb_files:
        try:
            structures[path] = load_structure(path, parser)
        except Exception as e:
            load_errors.append((path, str(e)))
            print(f"  [SKIP] {path}: {e}", file=sys.stderr)

    if not structures:
        print("Error: no structures loaded successfully.", file=sys.stderr)
        sys.exit(1)

    # Load reference if provided
    use_centroid = args.reference is None
    ref_struct = None
    if not use_centroid:
        try:
            ref_struct = load_structure(args.reference, parser)
            print(f"Reference: {args.reference}")
        except Exception as e:
            print(f"Error loading reference: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("No reference provided — will use centroid of all structures.")

    # Align all structures to reference (or first pass for centroid)
    # For centroid: align all to the first structure as a common frame,
    # then compute centroid coords in that frame.
    align_errors = {}

    if use_centroid:
        # Pick first successfully loaded structure as alignment anchor
        anchor_path = list(structures.keys())[0]
        anchor_struct = structures[anchor_path]
        print(f"Alignment anchor (for centroid): {anchor_path}")

        for path, struct in structures.items():
            if path == anchor_path:
                continue
            try:
                align_to_reference(struct, anchor_struct, args.chain)
            except Exception as e:
                align_errors[path] = str(e)
                print(f"  [ALIGN ERROR] {path}: {e}", file=sys.stderr)
    else:
        for path, struct in structures.items():
            try:
                align_to_reference(struct, ref_struct, args.chain)
            except Exception as e:
                align_errors[path] = str(e)
                print(f"  [ALIGN ERROR] {path}: {e}", file=sys.stderr)

    # Extract ligand coords
    ligand_coords = {}
    coord_errors = {}
    for path, struct in structures.items():
        if path in align_errors:
            continue
        try:
            ligand_coords[path] = get_ligand_heavy_atoms(struct, args.ligand, args.chain)
        except Exception as e:
            coord_errors[path] = str(e)
            print(f"  [LIGAND ERROR] {path}: {e}", file=sys.stderr)

    if not ligand_coords:
        print("Error: no ligand coordinates extracted.", file=sys.stderr)
        sys.exit(1)

    # Build reference coords
    if use_centroid:
        centroid = compute_centroid_coords(list(ligand_coords.values()))
        ref_coords = centroid
        n_centroid_atoms = len(centroid)
        print(f"Centroid built from {len(ligand_coords)} structures "
              f"over {n_centroid_atoms} common heavy atoms.")
    else:
        try:
            ref_coords = get_ligand_heavy_atoms(ref_struct, args.ligand, args.chain)
        except Exception as e:
            print(f"Error extracting ligand from reference: {e}", file=sys.stderr)
            sys.exit(1)

    # Compute RMSD for each structure
    rmsd_results = {}   # path -> (rmsd, n_atoms)
    rmsd_errors = {}
    for path, coords in ligand_coords.items():
        try:
            rmsd, n = compute_rmsd(coords, ref_coords)
            rmsd_results[path] = (rmsd, n)
        except Exception as e:
            rmsd_errors[path] = str(e)
            print(f"  [RMSD ERROR] {path}: {e}", file=sys.stderr)

    if not rmsd_results:
        print("Error: no RMSD values computed.", file=sys.stderr)
        sys.exit(1)

    # Z-score computation
    rmsd_values = np.array([v[0] for v in rmsd_results.values()])
    mean_rmsd = rmsd_values.mean()
    std_rmsd = rmsd_values.std()

    zscores = {}
    for path, (rmsd, _) in rmsd_results.items():
        zscores[path] = (rmsd - mean_rmsd) / std_rmsd if std_rmsd > 0 else 0.0

    # Determine inclusion
    included = {}
    for path in rmsd_results:
        z = zscores[path]
        if args.no_exclude:
            included[path] = True
        else:
            included[path] = abs(z) <= args.zscore

    # Final average over included structures
    included_rmsds = [rmsd_results[p][0] for p in rmsd_results if included[p]]
    final_avg = np.mean(included_rmsds) if included_rmsds else float("nan")
    n_included = len(included_rmsds)
    n_excluded = len(rmsd_results) - n_included

    # Build table
    rows = []
    for path in sorted(rmsd_results.keys()):
        rmsd, n_atoms = rmsd_results[path]
        z = zscores[path]
        inc = included[path]
        rows.append([
            Path(path).name,
            f"{rmsd:.3f}",
            f"{z:+.2f}",
            n_atoms,
            "Yes" if inc else "No (outlier)",
        ])

    # Add error rows
    all_error_paths = (
        set(load_errors) |
        set(align_errors) |
        set(coord_errors) |
        set(rmsd_errors)
    )
    for path, err in sorted(align_errors.items()):
        rows.append([Path(path).name, "—", "—", "—", f"Error: {err[:40]}"])
    for path, err in sorted(coord_errors.items()):
        rows.append([Path(path).name, "—", "—", "—", f"Error: {err[:40]}"])
    for path, err in sorted(rmsd_errors.items()):
        rows.append([Path(path).name, "—", "—", "—", f"Error: {err[:40]}"])

    headers = ["File", "RMSD (Å)", "Z-score", "Atoms matched", "Included"]

    print()
    print(tabulate(rows, headers=headers, tablefmt="github"))
    print()
    print(f"Ligand:            {args.ligand}")
    if not use_centroid:
        print(f"Reference:         {Path(args.reference).name}")
    else:
        print(f"Reference:         centroid ({n_centroid_atoms} atoms)")
    print(f"Z-score threshold: {args.zscore}")
    print(f"Structures:        {len(rmsd_results)} computed, "
          f"{n_included} included, {n_excluded} excluded")
    print(f"Average RMSD:      {final_avg:.3f} Å  (over {n_included} structures)")


if __name__ == "__main__":
    main()
