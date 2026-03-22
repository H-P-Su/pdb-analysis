#!/usr/bin/env python3
"""
visualize_interactions.py — Visualize structural interactions.

Generates two outputs per ligand:
  1. Interactive 3D HTML (py3Dmol) — rotatable, zoomable in any browser
  2. 2D interaction fingerprint PNG (matplotlib) — publication-ready summary

Color scheme:
  Ligand          — yellow carbons
  Contacts        — cyan sticks
  H-bonds         — orange sticks + yellow dashed lines
  Pi interactions — magenta sticks + magenta dashed lines
  Salt bridges    — blue/red sticks + orange dashed lines

All functions are importable for use in other projects.

CLI:
    python3 visualize_interactions.py structure.pdb --analyze "resn LIG"
    python3 visualize_interactions.py structure.pdb --analyze "organic" --outdir ./images
    python3 visualize_interactions.py structure.pdb --analyze "organic" --width 1200 --height 900
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from analyze_ligands import (
    Contact, HBond, LigandReport, PiInteraction, SaltBridge,
    SelectionParser, analyze_ligand, load_structure, primary_ligands,
    WATER_NAMES,
)
from Bio.PDB.Structure import Structure

# ─── Amino acid lookup ────────────────────────────────────────────────────────

_AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# ─── Palette ──────────────────────────────────────────────────────────────────

_COLORS = {
    "ligand":   {"scheme": "yellowCarbon",   "hex": "#e8c000"},
    "contact":  {"scheme": "cyanCarbon",     "hex": "#4dd4e4"},
    "hbond":    {"scheme": "orangeCarbon",   "hex": "#ff8c00"},
    "pi":       {"scheme": "magentaCarbon",  "hex": "#cc44cc"},
    "salt_cat": {"scheme": "blueCarbon",     "hex": "#4488ff"},
    "salt_ani": {"scheme": "redCarbon",      "hex": "#ff4444"},
}

_DASH_COLORS = {
    "hbond": "#ffdd00",
    "pi":    "#ff44ff",
    "salt":  "#ff8800",
}

_CPK_COLORS = {
    "C": "#404040", "N": "#3050F8", "O": "#FF0D0D", "S": "#FFFF30",
    "P": "#FF8000", "F": "#90E050", "CL": "#1FF01F", "BR": "#A62929",
    "I": "#940094", "H": "#FFFFFF", "FE": "#E06633", "ZN": "#7D80B0",
    "MG": "#8AFF00", "CA": "#3DFF00", "MN": "#9C7AC7", "CU": "#C88033",
}


def _draw_spoked_arc(ax, res_pos: np.ndarray, direction: np.ndarray,
                     outer_r: float = 0.45, n_spokes: int = 6) -> None:
    """Draw a LIGPLOT-style spoked semicircle indicating a hydrophobic contact."""
    angle_to_lig = np.degrees(np.arctan2(direction[1], direction[0]))
    arc_rad = np.radians(np.linspace(angle_to_lig - 70, angle_to_lig + 70, 60))
    ax.plot(res_pos[0] + outer_r * np.cos(arc_rad),
            res_pos[1] + outer_r * np.sin(arc_rad),
            "-", color="#999999", linewidth=1.4, zorder=1)
    inner_r = outer_r * 0.45
    for theta_deg in np.linspace(angle_to_lig - 60, angle_to_lig + 60, n_spokes):
        theta = np.radians(theta_deg)
        ax.plot([res_pos[0] + inner_r * np.cos(theta),
                 res_pos[0] + outer_r * np.cos(theta)],
                [res_pos[1] + inner_r * np.sin(theta),
                 res_pos[1] + outer_r * np.sin(theta)],
                "-", color="#999999", linewidth=0.8, zorder=1)


# ─── Geometry helpers ─────────────────────────────────────────────────────────

def _get_atom_coord(structure: Structure, chain: str,
                    resi: int, atom_name: str) -> np.ndarray | None:
    for model in structure:
        for ch in model:
            if ch.id != chain:
                continue
            for res in ch:
                if res.id[1] != resi:
                    continue
                if atom_name in res:
                    return res[atom_name].get_vector().get_array()
    return None


def _get_centroid(structure: Structure, chain: str,
                  resi: int, atom_names: list[str]) -> np.ndarray | None:
    coords = []
    for model in structure:
        for ch in model:
            if ch.id != chain:
                continue
            for res in ch:
                if res.id[1] != resi:
                    continue
                for name in atom_names:
                    if name in res:
                        coords.append(res[name].get_vector().get_array())
    return np.mean(coords, axis=0) if coords else None


def _add_dashed_line(view, p1: np.ndarray, p2: np.ndarray,
                     color: str, n_dashes: int = 7, radius: float = 0.07):
    """Simulate a dashed line with alternating cylinders (py3Dmol notebook use)."""
    p1, p2 = np.asarray(p1, float), np.asarray(p2, float)
    for i in range(n_dashes):
        t0 = (2 * i)     / (2 * n_dashes)
        t1 = (2 * i + 1) / (2 * n_dashes)
        t1 = min(t1, 1.0)
        s = (p1 + t0 * (p2 - p1)).tolist()
        e = (p1 + t1 * (p2 - p1)).tolist()
        view.addCylinder({
            "start": {"x": s[0], "y": s[1], "z": s[2]},
            "end":   {"x": e[0], "y": e[1], "z": e[2]},
            "radius": radius,
            "color": color,
            "fromCap": 1,
            "toCap": 1,
        })


def _dashed_line_specs(p1: np.ndarray, p2: np.ndarray,
                       color: str, n_dashes: int = 7,
                       radius: float = 0.07) -> list[dict]:
    """Return a list of 3Dmol.js cylinder spec dicts that form a dashed line."""
    p1, p2 = np.asarray(p1, float), np.asarray(p2, float)
    specs = []
    for i in range(n_dashes):
        t0 = (2 * i)     / (2 * n_dashes)
        t1 = min((2 * i + 1) / (2 * n_dashes), 1.0)
        s = (p1 + t0 * (p2 - p1)).tolist()
        e = (p1 + t1 * (p2 - p1)).tolist()
        specs.append({
            "start": {"x": s[0], "y": s[1], "z": s[2]},
            "end":   {"x": e[0], "y": e[1], "z": e[2]},
            "radius": radius, "color": color,
            "fromCap": 1, "toCap": 1,
        })
    return specs


# ─── 3D HTML viewer (py3Dmol) ─────────────────────────────────────────────────

def build_view(structure_path: str | Path,
               report: LigandReport,
               structure: Structure | None = None,
               width: int = 900,
               height: int = 700):
    """
    Build a py3Dmol view for a LigandReport.

    Parameters
    ----------
    structure_path :
        Path to the PDB or CIF file (read as text for py3Dmol).
    report :
        LigandReport from analyze_ligand().
    structure :
        Pre-loaded BioPython Structure.  Loaded from structure_path if None.
    width, height :
        Viewer dimensions in pixels.

    Returns
    -------
    Configured py3Dmol.view object — call .show() in a notebook or
    pass to save_html() to write an HTML file.
    """
    import py3Dmol
    structure_path = Path(structure_path)
    if structure is None:
        structure = load_structure(structure_path)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(structure_path.read_text(), "pdb")

    lig_key = (report.ligand_chain, report.ligand_resi)

    # ── Base: everything as gray cartoon ──────────────────────────────────────
    view.setStyle({}, {"cartoon": {"color": "#b8b8b8", "opacity": 0.75}})
    view.setStyle({"resn": list(WATER_NAMES)}, {})   # hide water

    # ── Ligand: sticks only (replaces cartoon) ────────────────────────────────
    view.setStyle(
        {"chain": report.ligand_chain, "resi": report.ligand_resi},
        {"stick": {"colorscheme": _COLORS["ligand"]["scheme"], "radius": 0.25}},
    )

    # ── Collect residue sets by interaction type ───────────────────────────────
    contact_res: set[tuple] = set()
    hbond_res:   set[tuple] = set()
    pi_res:      set[tuple] = set()
    sb_cat_res:  set[tuple] = set()
    sb_ani_res:  set[tuple] = set()

    for c in report.contacts:
        r = (c.chain2, c.resi2)
        if r != lig_key:
            contact_res.add(r)

    for h in report.hbonds:
        for chain, resi in [(h.donor_chain, h.donor_resi),
                             (h.acceptor_chain, h.acceptor_resi)]:
            if (chain, resi) != lig_key:
                hbond_res.add((chain, resi))

    for p in report.pi_interactions:
        for chain, resi in [(p.chain1, p.resi1), (p.chain2, p.resi2)]:
            if (chain, resi) != lig_key:
                pi_res.add((chain, resi))

    for s in report.salt_bridges:
        if (s.cation_chain, s.cation_resi) != lig_key:
            sb_cat_res.add((s.cation_chain, s.cation_resi))
        if (s.anion_chain, s.anion_resi) != lig_key:
            sb_ani_res.add((s.anion_chain, s.anion_resi))

    # ── Add sticks per type (higher priority overrides lower) ─────────────────
    def _add_sticks(res_set: set[tuple], scheme: str, radius: float = 0.18):
        for chain, resi in res_set:
            view.addStyle(
                {"chain": chain, "resi": resi},
                {"stick": {"colorscheme": scheme, "radius": radius}},
            )

    plain_contacts = contact_res - hbond_res - pi_res - sb_cat_res - sb_ani_res
    plain_pi       = pi_res - hbond_res
    plain_sb_cat   = sb_cat_res - hbond_res
    plain_sb_ani   = sb_ani_res - hbond_res

    _add_sticks(plain_contacts, _COLORS["contact"]["scheme"],  0.15)
    _add_sticks(plain_pi,       _COLORS["pi"]["scheme"],       0.20)
    _add_sticks(plain_sb_cat,   _COLORS["salt_cat"]["scheme"], 0.20)
    _add_sticks(plain_sb_ani,   _COLORS["salt_ani"]["scheme"], 0.20)
    _add_sticks(hbond_res,      _COLORS["hbond"]["scheme"],    0.20)  # highest priority

    # ── H-bond dashed lines ────────────────────────────────────────────────────
    for hb in report.hbonds:
        c1 = _get_atom_coord(structure, hb.donor_chain,    hb.donor_resi,    hb.donor_atom)
        c2 = _get_atom_coord(structure, hb.acceptor_chain, hb.acceptor_resi, hb.acceptor_atom)
        if c1 is not None and c2 is not None:
            _add_dashed_line(view, c1, c2, _DASH_COLORS["hbond"], radius=0.07)

    # ── Pi dashed lines (centroid-to-centroid) ─────────────────────────────────
    for pi in report.pi_interactions:
        atoms1 = pi.ring1_label.split(",")
        atoms2 = pi.ring2_label.split(",")
        if pi.subtype == "cation_pi":
            c1 = _get_atom_coord(structure, pi.chain1, pi.resi1, atoms1[0])
        else:
            c1 = _get_centroid(structure, pi.chain1, pi.resi1, atoms1)
        c2 = _get_centroid(structure, pi.chain2, pi.resi2, atoms2)
        if c1 is not None and c2 is not None:
            _add_dashed_line(view, c1, c2, _DASH_COLORS["pi"], radius=0.06)

    # ── Salt bridge dashed lines ───────────────────────────────────────────────
    for sb in report.salt_bridges:
        c1 = _get_atom_coord(structure, sb.cation_chain, sb.cation_resi, sb.cation_atom)
        c2 = _get_atom_coord(structure, sb.anion_chain,  sb.anion_resi,  sb.anion_atom)
        if c1 is not None and c2 is not None:
            _add_dashed_line(view, c1, c2, _DASH_COLORS["salt"], radius=0.06)

    # ── Labels for H-bond, pi, and salt bridge residues ───────────────────────
    labeled = hbond_res | pi_res | sb_cat_res | sb_ani_res
    for chain, resi in labeled:
        # Try CA first; fall back to first heavy atom in the residue
        coord = _get_atom_coord(structure, chain, resi, "CA")
        if coord is None:
            coord = _get_atom_coord(structure, chain, resi, "CB")
        if coord is None:
            continue
        resn = _get_resname(structure, chain, resi)
        color = ("#ff8c00" if (chain, resi) in hbond_res
                 else "#cc44cc" if (chain, resi) in pi_res
                 else "#4488ff" if (chain, resi) in sb_cat_res
                 else "#ff4444")
        view.addLabel(
            f"{resn}{resi}",
            {"backgroundColor": color, "fontColor": "white",
             "fontSize": 11, "showBackground": True, "borderThickness": 0},
            {"chain": chain, "resi": resi, "atom": "CA"},
        )

    view.setBackgroundColor("white")
    view.zoomTo({"chain": report.ligand_chain, "resi": report.ligand_resi})

    return view


def _get_resname(structure: Structure, chain: str, resi: int) -> str:
    for model in structure:
        for ch in model:
            if ch.id != chain:
                continue
            for res in ch:
                if res.id[1] == resi:
                    return res.get_resname().strip()
    return "?"


def save_html(structure_path: str | Path,
              report: LigandReport,
              output_path: str | Path,
              structure: Structure | None = None,
              width: int = 900,
              height: int = 700,
              active_types: set[str] | None = None) -> Path:
    """
    Save an interactive 3D HTML viewer for a LigandReport.

    The HTML file is self-contained (loads 3Dmol.js via CDN) and can be
    opened in any modern browser.

    Parameters
    ----------
    structure_path :
        Path to the structure file.
    report :
        LigandReport from analyze_ligand().
    output_path :
        Output .html file path.
    structure :
        Pre-loaded BioPython Structure (optional).
    width, height :
        Viewer dimensions in pixels.

    Returns
    -------
    Path to the saved HTML file.
    """
    structure_path = Path(structure_path)
    if structure is None:
        structure = load_structure(structure_path)
    html = _render_html(structure_path, report, structure, width, height, active_types)
    output_path = Path(output_path)
    output_path.write_text(html)
    print(f"  Saved 3D viewer    → {output_path}")
    return output_path


def _render_html(structure_path: Path, report: LigandReport,
                 structure: Structure, width: int, height: int,
                 active_types: set[str] | None = None) -> str:
    """
    Generate a self-contained HTML page with:
      - 3Dmol.js interactive viewer (loads from CDN, no py3Dmol needed)
      - Per-type toggle buttons to show/hide contacts, H-bonds, pi, salt bridges
      - Residue panel listing every interaction partner with atom names and distances
    """
    import json

    lig_key = (report.ligand_chain, report.ligand_resi)
    pdb_text = Path(structure_path).read_text()

    # ── Per-type atom style lists  {sel, style} ──────────────────────────────
    # These are applied in priority order in JS: contacts < pi < salt < hbonds
    contact_atom_styles: list[dict] = []
    hbond_atom_styles:   list[dict] = []
    pi_atom_styles:      list[dict] = []
    salt_atom_styles:    list[dict] = []

    # ── Per-type shape specs (dashed lines) ───────────────────────────────────
    hbond_shape_specs: list[dict] = []
    pi_shape_specs:    list[dict] = []
    salt_shape_specs:  list[dict] = []

    # ── Per-type label specs ──────────────────────────────────────────────────
    hbond_label_specs: list[dict] = []
    pi_label_specs:    list[dict] = []
    salt_label_specs:  list[dict] = []

    # ── Side-panel HTML lists ─────────────────────────────────────────────────
    contact_items: list[str] = []
    hbond_items:   list[str] = []
    pi_items:      list[str] = []
    salt_items:    list[str] = []

    # helper — deduplicated atom style append
    def _atom_style(bucket: list, chain: str, resi: int, scheme: str, radius: float):
        bucket.append({"sel": {"chain": chain, "resi": resi},
                       "style": {"stick": {"colorscheme": scheme, "radius": radius}}})

    def _label_spec(chain: str, resi: int, resn: str, bg: str) -> Optional[dict]:
        coord = _get_atom_coord(structure, chain, resi, "CA")
        if coord is None:
            coord = _get_atom_coord(structure, chain, resi, "CB")
        if coord is None:
            return None
        return {"text": f"{resn}{resi}",
                "style": {"backgroundColor": bg, "fontColor": "white",
                          "fontSize": 11, "showBackground": True,
                          "borderThickness": 0},
                "sel": {"chain": chain, "resi": resi, "atom": "CA"}}

    # ── Contacts ──────────────────────────────────────────────────────────────
    contact_residues: dict[tuple, dict] = {}
    for c in report.contacts:
        key = (c.chain2, c.resi2)
        if key == lig_key:
            continue
        if key not in contact_residues:
            contact_residues[key] = {"resn": c.resn2, "dists": []}
        contact_residues[key]["dists"].append(c.distance)

    seen_contact: set[tuple] = set()
    for (chain, resi), info in sorted(contact_residues.items(), key=lambda x: x[0][1]):
        if (chain, resi) not in seen_contact:
            _atom_style(contact_atom_styles, chain, resi,
                        _COLORS["contact"]["scheme"], 0.15)
            seen_contact.add((chain, resi))
        n = len(info["dists"])
        min_d = min(info["dists"])
        contact_items.append(
            f'<li><span class="rn">{info["resn"]}{resi}:{chain}</span>'
            f'<span class="dt">{n} contact{"s" if n>1 else ""} &nbsp;·&nbsp; min {min_d:.1f}&thinsp;Å</span></li>'
        )

    # ── H-bonds ───────────────────────────────────────────────────────────────
    seen_hbond: set[tuple] = set()
    hbond_label_seen: set[tuple] = set()
    for hb in report.hbonds:
        c1 = _get_atom_coord(structure, hb.donor_chain,    hb.donor_resi,    hb.donor_atom)
        c2 = _get_atom_coord(structure, hb.acceptor_chain, hb.acceptor_resi, hb.acceptor_atom)
        if c1 is not None and c2 is not None:
            hbond_shape_specs.extend(
                _dashed_line_specs(c1, c2, _DASH_COLORS["hbond"], radius=0.07))

        for chain, resi, resn, scheme in [
            (hb.donor_chain,    hb.donor_resi,    hb.donor_resn,
             _COLORS["hbond"]["scheme"]),
            (hb.acceptor_chain, hb.acceptor_resi, hb.acceptor_resn,
             _COLORS["hbond"]["scheme"]),
        ]:
            key = (chain, resi)
            if key == lig_key:
                continue
            if key not in seen_hbond:
                _atom_style(hbond_atom_styles, chain, resi, scheme, 0.20)
                seen_hbond.add(key)
            if key not in hbond_label_seen:
                spec = _label_spec(chain, resi, resn, _COLORS["hbond"]["hex"])
                if spec:
                    hbond_label_specs.append(spec)
                hbond_label_seen.add(key)

        # Build display strings for donor and acceptor sides
        d_lig = (hb.donor_chain    == report.ligand_chain and
                 hb.donor_resi     == report.ligand_resi)
        a_lig = (hb.acceptor_chain == report.ligand_chain and
                 hb.acceptor_resi  == report.ligand_resi)
        if d_lig:
            donor_str    = f'<span class="lig">{report.ligand_resn} ({hb.donor_atom})</span>'
            acceptor_str = (f'<span class="rn">{hb.acceptor_resn}{hb.acceptor_resi}'
                            f':{hb.acceptor_chain}</span>'
                            f' <span class="at">({hb.acceptor_atom})</span>')
        elif a_lig:
            donor_str    = (f'<span class="rn">{hb.donor_resn}{hb.donor_resi}'
                            f':{hb.donor_chain}</span>'
                            f' <span class="at">({hb.donor_atom})</span>')
            acceptor_str = f'<span class="lig">{report.ligand_resn} ({hb.acceptor_atom})</span>'
        else:
            donor_str    = (f'<span class="rn">{hb.donor_resn}{hb.donor_resi}'
                            f':{hb.donor_chain}</span>'
                            f' <span class="at">({hb.donor_atom})</span>')
            acceptor_str = (f'<span class="rn">{hb.acceptor_resn}{hb.acceptor_resi}'
                            f':{hb.acceptor_chain}</span>'
                            f' <span class="at">({hb.acceptor_atom})</span>')
        hbond_items.append(
            f'<li>{donor_str} <span class="arrow">→</span> {acceptor_str}'
            f' <span class="dt">{hb.distance:.1f}&thinsp;Å</span></li>'
        )

    # ── Pi interactions ───────────────────────────────────────────────────────
    seen_pi: set[tuple] = set()
    pi_label_seen: set[tuple] = set()
    for pi in report.pi_interactions:
        atoms1 = pi.ring1_label.split(",")
        atoms2 = pi.ring2_label.split(",")
        c1 = (_get_atom_coord(structure, pi.chain1, pi.resi1, atoms1[0])
              if pi.subtype == "cation_pi"
              else _get_centroid(structure, pi.chain1, pi.resi1, atoms1))
        c2 = _get_centroid(structure, pi.chain2, pi.resi2, atoms2)
        if c1 is not None and c2 is not None:
            pi_shape_specs.extend(
                _dashed_line_specs(c1, c2, _DASH_COLORS["pi"], radius=0.06))

        for chain, resi, resn in [(pi.chain1, pi.resi1, pi.resn1),
                                   (pi.chain2, pi.resi2, pi.resn2)]:
            key = (chain, resi)
            if key == lig_key:
                continue
            if key not in seen_pi:
                _atom_style(pi_atom_styles, chain, resi,
                            _COLORS["pi"]["scheme"], 0.20)
                seen_pi.add(key)
            if key not in pi_label_seen:
                spec = _label_spec(chain, resi, resn, _COLORS["pi"]["hex"])
                if spec:
                    pi_label_specs.append(spec)
                pi_label_seen.add(key)

        subtype_label = pi.subtype.replace("_", "-")
        # Identify the protein residue (non-ligand side)
        if (pi.chain1, pi.resi1) != lig_key:
            prot_str = (f'<span class="rn">{pi.resn1}{pi.resi1}:{pi.chain1}</span>')
        else:
            prot_str = (f'<span class="rn">{pi.resn2}{pi.resi2}:{pi.chain2}</span>')
        pi_items.append(
            f'<li>{prot_str} <span class="at">{subtype_label}</span>'
            f' <span class="dt">{pi.center_distance:.1f}&thinsp;Å</span></li>'
        )

    # ── Salt bridges ──────────────────────────────────────────────────────────
    seen_salt: set[tuple] = set()
    salt_label_seen: set[tuple] = set()
    for sb in report.salt_bridges:
        c1 = _get_atom_coord(structure, sb.cation_chain, sb.cation_resi, sb.cation_atom)
        c2 = _get_atom_coord(structure, sb.anion_chain,  sb.anion_resi,  sb.anion_atom)
        if c1 is not None and c2 is not None:
            salt_shape_specs.extend(
                _dashed_line_specs(c1, c2, _DASH_COLORS["salt"], radius=0.06))

        for chain, resi, resn, scheme, bg, role in [
            (sb.cation_chain, sb.cation_resi, sb.cation_resn,
             _COLORS["salt_cat"]["scheme"], _COLORS["salt_cat"]["hex"], "+"),
            (sb.anion_chain,  sb.anion_resi,  sb.anion_resn,
             _COLORS["salt_ani"]["scheme"], _COLORS["salt_ani"]["hex"], "−"),
        ]:
            key = (chain, resi)
            if key == lig_key:
                continue
            if key not in seen_salt:
                _atom_style(salt_atom_styles, chain, resi, scheme, 0.20)
                seen_salt.add(key)
            if key not in salt_label_seen:
                spec = _label_spec(chain, resi, resn, bg)
                if spec:
                    salt_label_specs.append(spec)
                salt_label_seen.add(key)

        salt_items.append(
            f'<li>'
            f'<span class="rn">{sb.cation_resn}{sb.cation_resi}:{sb.cation_chain}</span>'
            f' <span class="at role-pos">(+) {sb.cation_atom}</span>'
            f' <span class="arrow">↔</span> '
            f'<span class="at role-neg">{sb.anion_atom} (−)</span>'
            f'<span class="rn">{sb.anion_resn}{sb.anion_resi}:{sb.anion_chain}</span>'
            f' <span class="dt">{sb.distance:.1f}&thinsp;Å</span></li>'
        )

    # ── Side panel HTML sections ──────────────────────────────────────────────
    def _section(title: str, color: str, items: list[str], tid: str) -> str:
        count = f" ({len(items)})" if items else ""
        body  = (f"<ul>{''.join(items)}</ul>" if items
                 else '<p class="empty">None detected</p>')
        return (f'<div class="isect" id="list-{tid}">'
                f'<div class="isect-title" style="color:{color}">{title}{count}</div>'
                f'{body}</div>')

    panel_html = (
        _section("Contacts",       _COLORS["contact"]["hex"],  contact_items, "contacts") +
        _section("H-bonds",        _COLORS["hbond"]["hex"],    hbond_items,   "hbonds")   +
        _section("Pi interactions",_COLORS["pi"]["hex"],       pi_items,      "pi")       +
        _section("Salt bridges",   _COLORS["salt_cat"]["hex"], salt_items,    "salt")
    )

    # ── Sequence data ─────────────────────────────────────────────────────────
    seq_residues: list[dict] = []
    seen_seq_keys: set[tuple] = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] != ' ':  # skip HETATM and water
                    continue
                resn  = residue.get_resname().strip()
                if resn not in _AA3TO1:
                    continue
                resi  = residue.get_id()[1]
                cid   = chain.get_id()
                key   = (cid, resi)
                if key not in seen_seq_keys:
                    seen_seq_keys.add(key)
                    seq_residues.append({"resi": resi, "resn": resn,
                                         "aa": _AA3TO1[resn], "chain": cid})
        break  # first model only
    seq_js = json.dumps(seq_residues)

    # ── Header subtitle ───────────────────────────────────────────────────────
    subtitle = (f"{report.ligand_resn} {report.ligand_chain}:{report.ligand_resi}"
                f" &nbsp;·&nbsp; {report.n_contacts} contacts"
                f" &nbsp;·&nbsp; {report.n_hbonds} H-bonds"
                f" &nbsp;·&nbsp; {report.n_pi} pi"
                f" &nbsp;·&nbsp; {report.n_salt_bridges} salt bridges")

    # ── Serialise data for JS ─────────────────────────────────────────────────
    J = json.dumps
    water_js = J(sorted(WATER_NAMES))
    # active_types controls which checkboxes start checked (None = all checked)
    _all = active_types is None
    def _chk(t: str) -> str:
        return "checked" if (_all or t in active_types) else ""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{report.ligand_resn} Interactions</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
        background: #f0f2f5; color: #222; height: 100vh; display: flex;
        flex-direction: column; overflow: hidden; }}
.header {{ background: #1e2d3e; color: white; padding: 10px 18px;
           flex-shrink: 0; }}
.header h2 {{ font-size: 15px; margin-bottom: 3px; }}
.subtitle  {{ font-size: 12px; opacity: 0.8; }}
.body {{ display: flex; flex: 1; overflow: hidden; }}
#viewer-wrap {{ flex: 1; position: relative; background: white; }}
#viewer {{ width: 100%; height: 100%; position: absolute; }}
.panel {{ width: 290px; flex-shrink: 0; background: white;
          border-left: 1px solid #dde; display: flex; flex-direction: column;
          overflow: hidden; }}
.toggles {{ padding: 10px 12px; border-bottom: 1px solid #eee; flex-shrink: 0; }}
.toggles-title {{ font-size: 11px; font-weight: 700; letter-spacing: .06em;
                  color: #888; text-transform: uppercase; margin-bottom: 8px; }}
.tog-row {{ display: flex; align-items: center; gap: 8px;
            margin-bottom: 5px; cursor: pointer; }}
.tog-row input {{ cursor: pointer; accent-color: #4a90d9; width: 15px; height: 15px; }}
.tog-label {{ font-size: 13px; display: flex; align-items: center; gap: 5px; }}
.tog-dot {{ width: 11px; height: 11px; border-radius: 50%; flex-shrink: 0; }}
.tog-actions {{ display: flex; gap: 6px; margin-top: 8px; }}
.tog-btn {{ font-size: 11px; padding: 2px 8px; border: 1px solid #ccc;
            border-radius: 3px; cursor: pointer; background: #f7f7f7; }}
.tog-btn:hover {{ background: #eee; }}
.residues {{ flex: 1; overflow-y: auto; padding: 8px 12px; }}
.isect {{ margin-bottom: 12px; }}
.isect-title {{ font-size: 12px; font-weight: 700; margin-bottom: 5px;
                padding-bottom: 3px; border-bottom: 2px solid currentColor; }}
.isect ul {{ list-style: none; padding: 0; }}
.isect li {{ font-size: 12px; padding: 3px 0; border-bottom: 1px solid #f2f2f2;
             line-height: 1.5; }}
.rn   {{ font-weight: 600; font-family: monospace; }}
.lig  {{ font-weight: 600; font-family: monospace; color: #b89000; }}
.at   {{ color: #666; font-size: 11px; }}
.dt   {{ color: #888; font-size: 11px; float: right; }}
.arrow {{ color: #aaa; }}
.role-pos {{ color: #4488ff; }}
.role-neg {{ color: #ff4444; }}
.empty {{ font-size: 12px; color: #aaa; font-style: italic; }}
#vi-seq-panel {{
  background:#1c2128; border-bottom:1px solid #30363d;
  padding:6px 12px; white-space:pre-wrap; word-break:break-all;
  font-family:monospace; font-size:12px; line-height:1.7; user-select:none;
  flex-shrink:0; max-height:90px; overflow-y:auto;
}}
.vi-seq-aa {{ cursor:pointer; padding:0 1px; border-radius:2px; transition:background 0.1s; color:#cdd3de; }}
.vi-seq-aa:hover {{ background:#30363d; }}
.vi-seq-sep {{ color:#8b949e; cursor:default; }}
.vi-seq-chain-lbl {{ color:#8b949e; font-size:10px; font-family:sans-serif; font-weight:600;
                    cursor:default; letter-spacing:0.5px; }}
#vi-controls {{
  position:absolute; top:10px; left:10px; z-index:100;
  background:rgba(13,17,23,0.88); border-radius:8px; padding:10px 14px;
  font-size:13px; min-width:220px; color:#cdd3de;
}}
#vi-controls h4 {{ margin:0; color:#e6edf3; font-size:13px; }}
.vi-ctrl-header {{ display:flex; justify-content:space-between; align-items:center;
                  margin-bottom:8px; }}
.vi-collapse {{ background:none; border:none; color:#8b949e; cursor:pointer;
               font-size:13px; padding:0 2px; line-height:1; }}
.vi-collapse:hover {{ color:#e6edf3; }}
.vi-section {{ color:#8b949e; font-size:10px; text-transform:uppercase;
              letter-spacing:0.8px; margin:8px 0 4px; border-top:1px solid #21262d;
              padding-top:6px; }}
.vi-row {{ margin:3px 0; display:flex; align-items:center; gap:6px; cursor:pointer; }}
.vi-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:2px 6px; margin-bottom:2px; }}
.vi-bg-row {{ display:flex; align-items:center; gap:5px; margin-top:5px; }}
.vi-bg-btn {{ width:18px; height:18px; border-radius:3px; border:2px solid transparent;
             cursor:pointer; flex-shrink:0; }}
.vi-bg-btn.active {{ border-color:#e6edf3; }}
.vi-btn {{ background:rgba(255,255,255,0.08); border:1px solid #30363d; color:#cdd3de;
          border-radius:4px; padding:3px 9px; cursor:pointer; font-size:12px; }}
.vi-btn:hover {{ background:rgba(255,255,255,0.18); }}
#vi-toast {{ display:none; position:absolute; bottom:14px; right:14px; z-index:200;
            background:rgba(46,160,67,0.9); color:white; border-radius:6px;
            padding:5px 12px; font-size:12px; pointer-events:none; }}
</style>
</head>
<body>
<div class="header">
  <h2>Interaction Viewer</h2>
  <div class="subtitle">{subtitle}</div>
</div>
<div id="vi-seq-panel"><span id="vi-seq-letters"></span></div>
<div class="body">
  <div id="viewer-wrap"><div id="viewer"></div>
<div id="vi-controls">
  <div class="vi-ctrl-header">
    <h4>🔬 Viewer Controls</h4>
    <button class="vi-collapse" onclick="viToggleControls()" title="Collapse/expand" id="vi-toggle">▲</button>
  </div>
  <div id="vi-ctrl-body">
  <div class="vi-section">Style</div>
  <div class="vi-grid">
    <label class="vi-row"><input type="radio" name="vi-style" value="cartoon" checked><span>Cartoon</span></label>
    <label class="vi-row"><input type="radio" name="vi-style" value="stick"><span>Stick</span></label>
    <label class="vi-row"><input type="radio" name="vi-style" value="sphere"><span>Sphere</span></label>
    <label class="vi-row"><input type="radio" name="vi-style" value="line"><span>Line</span></label>
    <label class="vi-row"><input type="radio" name="vi-style" value="surface"><span>Surface</span></label>
  </div>
  <div id="vi-surf-opacity-row" style="display:none;margin:4px 0 2px;align-items:center;gap:6px;font-size:12px;color:#8b949e;">
    <span>Opacity</span>
    <input type="range" id="vi-surf-opacity" min="0" max="1" step="0.05" value="0.85"
           style="flex:1;accent-color:#58a6ff;" oninput="updateScene()">
    <span id="vi-surf-opacity-val">0.85</span>
  </div>

  <div class="vi-section">Color</div>
  <div class="vi-grid">
    <label class="vi-row"><input type="radio" name="vi-color" value="default" checked><span>Default</span></label>
    <label class="vi-row"><input type="radio" name="vi-color" value="spectrum"><span>Spectrum</span></label>
    <label class="vi-row"><input type="radio" name="vi-color" value="ss"><span>Sec. struct.</span></label>
    <label class="vi-row"><input type="radio" name="vi-color" value="chain"><span>Chain</span></label>
  </div>

  <div class="vi-section">Options</div>
  <label class="vi-row"><input type="checkbox" id="vi-spin"><span>Spin</span></label>
  <label class="vi-row"><input type="checkbox" id="vi-ortho"><span>Orthographic</span></label>
  <div class="vi-bg-row">
    <span style="color:#8b949e;font-size:11px;">BG</span>
    <div class="vi-bg-btn active" id="vi-bg0" style="background:#ffffff;border-color:#999;" title="White" onclick="viBg('#ffffff','vi-bg0')"></div>
    <div class="vi-bg-btn" id="vi-bg1" style="background:#0d1117;" title="Dark" onclick="viBg('#0d1117','vi-bg1')"></div>
    <div class="vi-bg-btn" id="vi-bg2" style="background:#000000;" title="Black" onclick="viBg('#000000','vi-bg2')"></div>
    <div class="vi-bg-btn" id="vi-bg3" style="background:#3d444d;" title="Grey" onclick="viBg('#3d444d','vi-bg3')"></div>
  </div>

  <div class="vi-section">Screenshot</div>
  <label class="vi-row"><input type="checkbox" id="vi-transparent"><span>Transparent BG</span></label>
  <div style="display:flex;gap:5px;margin-top:5px;">
    <button class="vi-btn" onclick="viSave()" title="Download PNG">📥 Save</button>
    <button class="vi-btn" onclick="viCopy()" title="Copy to clipboard">📋 Copy</button>
  </div>
  </div><!-- /vi-ctrl-body -->
</div>
<div id="vi-toast">✓ Copied to clipboard</div>
</div>
  <div class="panel">
    <div class="toggles">
      <div class="toggles-title">Show / Hide</div>
      <label class="tog-row">
        <input type="checkbox" id="tog-contacts" {_chk("contacts")} onchange="updateScene()">
        <span class="tog-label">
          <span class="tog-dot" style="background:{_COLORS['contact']['hex']}"></span>
          Contacts ({len(contact_items)})
        </span>
      </label>
      <label class="tog-row">
        <input type="checkbox" id="tog-hbonds" {_chk("hbonds")} onchange="updateScene()">
        <span class="tog-label">
          <span class="tog-dot" style="background:{_COLORS['hbond']['hex']}"></span>
          H-bonds ({len(hbond_items)})
        </span>
      </label>
      <label class="tog-row">
        <input type="checkbox" id="tog-pi" {_chk("pi")} onchange="updateScene()">
        <span class="tog-label">
          <span class="tog-dot" style="background:{_COLORS['pi']['hex']}"></span>
          Pi interactions ({len(pi_items)})
        </span>
      </label>
      <label class="tog-row">
        <input type="checkbox" id="tog-salt" {_chk("salt")} onchange="updateScene()">
        <span class="tog-label">
          <span class="tog-dot" style="background:{_COLORS['salt_cat']['hex']}"></span>
          Salt bridges ({len(salt_items)})
        </span>
      </label>
      <div class="tog-actions">
        <button class="tog-btn" onclick="setAll(true)">Show all</button>
        <button class="tog-btn" onclick="setAll(false)">Hide all</button>
      </div>
    </div>
    <div class="residues">{panel_html}</div>
  </div>
</div>
<script>
var PDB_DATA   = {J(pdb_text)};
var LIG_CHAIN  = {J(report.ligand_chain)};
var LIG_RESI   = {report.ligand_resi};
var WATER_NAMES= {water_js};

var atomStyles = {{
  contacts: {J(contact_atom_styles)},
  hbonds:   {J(hbond_atom_styles)},
  pi:       {J(pi_atom_styles)},
  salt:     {J(salt_atom_styles)}
}};
var shapeSpecs = {{
  hbonds: {J(hbond_shape_specs)},
  pi:     {J(pi_shape_specs)},
  salt:   {J(salt_shape_specs)}
}};
var labelSpecs = {{
  hbonds: {J(hbond_label_specs)},
  pi:     {J(pi_label_specs)},
  salt:   {J(salt_label_specs)}
}};

var viCurrentSurface = null;
var viSurfaceEpoch   = 0;
var viCurrentBg      = 'white';
var viSeqData        = {seq_js};
var viSelectedResi   = null;
var viSelectedChain  = null;

var VI_INT_COLORS = {{
  contacts: '{_COLORS["contact"]["hex"]}',
  hbonds:   '{_COLORS["hbond"]["hex"]}',
  pi:       '{_COLORS["pi"]["hex"]}',
  salt:     '{_COLORS["salt_cat"]["hex"]}'
}};

var viewer;

function viColorSeq() {{
  // Pre-compute all interaction types for each residue key "chain:resi"
  var residueTypes = {{}};
  var ALL_TYPES = ['contacts','pi','salt','hbonds'];
  ALL_TYPES.forEach(function(t) {{
    atomStyles[t].forEach(function(item) {{
      var key = item.sel.chain + ':' + item.sel.resi;
      if (!residueTypes[key]) residueTypes[key] = [];
      residueTypes[key].push(t);
    }});
  }});

  // Priority order (increasing — hbonds wins over salt wins over pi wins over contacts)
  var priority = ['contacts','pi','salt','hbonds'];

  document.querySelectorAll('.vi-seq-aa').forEach(function(el) {{
    var key = el.dataset.chain + ':' + el.dataset.resi;
    var resi = parseInt(el.dataset.resi);
    if (resi === viSelectedResi && el.dataset.chain === viSelectedChain) {{
      el.style.background = '#4488ff'; el.style.color = 'white'; return;
    }}
    var types = residueTypes[key] || [];
    var color = null;
    priority.forEach(function(t) {{
      var chk = document.getElementById('tog-' + t);
      if (chk && chk.checked && types.indexOf(t) >= 0) color = VI_INT_COLORS[t];
    }});
    if (color) {{
      el.style.background = color + '33';
      el.style.color = color;
    }} else {{
      el.style.background = ''; el.style.color = '';
    }}
  }});
}}

function viGetColorSpec() {{
  var s = document.querySelector('input[name="vi-color"]:checked').value;
  if (s === 'default')   return {{color: '#b8b8b8', opacity: 0.75}};
  if (s === 'spectrum')  return {{colorscheme: 'spectrum'}};
  if (s === 'ss')        return {{colorscheme: 'ssJmol'}};
  if (s === 'chain')     return {{colorscheme: 'chainHetatm'}};
  return {{color: '#b8b8b8', opacity: 0.75}};
}}

function viToggleControls() {{
  var body = document.getElementById('vi-ctrl-body');
  var btn = document.getElementById('vi-toggle');
  var collapsed = body.style.display === 'none';
  body.style.display = collapsed ? '' : 'none';
  btn.textContent = collapsed ? '▲' : '▼';
}}

function viBg(color, btnId) {{
  viCurrentBg = color;
  viewer.setBackgroundColor(color);
  viewer.render();
  document.querySelectorAll('.vi-bg-btn').forEach(function(b) {{ b.classList.remove('active'); }});
  document.getElementById(btnId).classList.add('active');
}}

function viParseBgRGB(hex) {{
  var c = document.createElement('canvas'); c.width = c.height = 1;
  var ctx = c.getContext('2d');
  ctx.fillStyle = hex; ctx.fillRect(0,0,1,1);
  return ctx.getImageData(0,0,1,1).data;
}}

function viCaptureURI() {{
  var baseURI = viewer.pngURI();
  if (!document.getElementById('vi-transparent').checked) return Promise.resolve(baseURI);
  return new Promise(function(resolve) {{
    var img = new Image();
    img.onload = function() {{
      var c = document.createElement('canvas');
      c.width = img.width; c.height = img.height;
      var ctx = c.getContext('2d');
      ctx.drawImage(img, 0, 0);
      var d = ctx.getImageData(0, 0, c.width, c.height);
      var px = d.data;
      var bg = viParseBgRGB(viCurrentBg);
      var thr = 18;
      for (var i = 0; i < px.length; i += 4) {{
        if (Math.abs(px[i]-bg[0]) + Math.abs(px[i+1]-bg[1]) + Math.abs(px[i+2]-bg[2]) < thr*3)
          px[i+3] = 0;
      }}
      ctx.putImageData(d, 0, 0);
      resolve(c.toDataURL('image/png'));
    }};
    img.src = baseURI;
  }});
}}

function viSave() {{
  viCaptureURI().then(function(uri) {{
    var a = document.createElement('a');
    a.href = uri; a.download = 'structure.png'; a.click();
  }});
}}

function viCopy() {{
  viCaptureURI().then(function(uri) {{
    fetch(uri).then(function(r) {{ return r.blob(); }}).then(function(blob) {{
      navigator.clipboard.write([new ClipboardItem({{'image/png': blob}})])
        .then(function() {{
          var t = document.getElementById('vi-toast');
          t.style.display = 'block';
          setTimeout(function() {{ t.style.display = 'none'; }}, 1800);
        }})
        .catch(function() {{
          alert('Clipboard blocked — use Save, or try in Chrome/Edge over HTTPS.');
        }});
    }});
  }});
}}

function updateScene() {{
  // 1. Reset: base style (respects vi-style and vi-color controls)
  if (viCurrentSurface !== null) {{
    viewer.removeSurface(viCurrentSurface); viCurrentSurface = null;
  }}
  viSurfaceEpoch++;
  var myEpoch = viSurfaceEpoch;

  var viStyleVal = document.querySelector('input[name="vi-style"]:checked').value;
  var viCol = viGetColorSpec();
  var viOpacityRow = document.getElementById('vi-surf-opacity-row');
  viOpacityRow.style.display = (viStyleVal === 'surface') ? 'flex' : 'none';

  viewer.setStyle({{}}, {{}});
  viewer.setStyle({{resn: WATER_NAMES}}, {{}});

  if (viStyleVal === 'cartoon') {{
    viewer.setStyle({{}}, {{cartoon: viCol}});
  }} else if (viStyleVal === 'stick') {{
    viewer.setStyle({{}}, {{stick: Object.assign({{radius: 0.12}}, viCol)}});
  }} else if (viStyleVal === 'sphere') {{
    viewer.setStyle({{}}, {{sphere: Object.assign({{scale: 0.4}}, viCol)}});
  }} else if (viStyleVal === 'line') {{
    viewer.setStyle({{}}, {{line: viCol}});
  }} else if (viStyleVal === 'surface') {{
    var viSurfOpacity = parseFloat(document.getElementById('vi-surf-opacity').value);
    document.getElementById('vi-surf-opacity-val').textContent = viSurfOpacity.toFixed(2);
    viewer.setStyle({{}}, {{cartoon: {{opacity: 0.2, color: '#888888'}}}});
    var surfP = viewer.addSurface($3Dmol.SurfaceType.VDW, Object.assign({{}}, viCol, {{opacity: viSurfOpacity}}), {{hetflag: false}});
    Promise.resolve(surfP).then(function(id) {{
      if (myEpoch !== viSurfaceEpoch) {{
        viewer.removeSurface(id);
      }} else {{
        viCurrentSurface = id;
      }}
    }});
  }}
  viewer.setStyle({{resn: WATER_NAMES}}, {{}});

  // 2. Ligand always on
  viewer.setStyle({{chain: LIG_CHAIN, resi: LIG_RESI}},
                  {{stick: {{colorscheme: "yellowCarbon", radius: 0.25}}}});

  // 3. Atom styles in ascending priority (hbonds last → wins)
  ["contacts","pi","salt","hbonds"].forEach(function(t) {{
    if (document.getElementById("tog-"+t).checked) {{
      atomStyles[t].forEach(function(item) {{
        viewer.addStyle(item.sel, item.style);
      }});
    }}
  }});

  // 4. Shapes (dashed lines)
  viewer.removeAllShapes();
  ["hbonds","pi","salt"].forEach(function(t) {{
    if (document.getElementById("tog-"+t).checked) {{
      shapeSpecs[t].forEach(function(spec) {{
        viewer.addCylinder(spec);
      }});
    }}
  }});

  // 5. Labels
  viewer.removeAllLabels();
  ["hbonds","pi","salt"].forEach(function(t) {{
    if (document.getElementById("tog-"+t).checked) {{
      labelSpecs[t].forEach(function(spec) {{
        viewer.addLabel(spec.text, spec.style, spec.sel);
      }});
    }}
  }});

  viewer.render();

  // 6. Side panel sections
  ["contacts","hbonds","pi","salt"].forEach(function(t) {{
    var el = document.getElementById("list-"+t);
    if (el) el.style.display =
      document.getElementById("tog-"+t).checked ? "" : "none";
  }});

  // 7. Sequence coloring
  viColorSeq();
}}

function setAll(on) {{
  ["contacts","hbonds","pi","salt"].forEach(function(t) {{
    document.getElementById("tog-"+t).checked = on;
  }});
  updateScene();
}}

// ── Build sequence strip ────────────────────────────────────────────────────
(function() {{
  var multiChain = viSeqData.length > 1 && viSeqData.some(function(r) {{
    return r.chain !== viSeqData[0].chain;
  }});
  var html = '';
  var prevChain = null;
  var chainPos  = 0;
  viSeqData.forEach(function(r) {{
    if (r.chain !== prevChain) {{
      if (prevChain !== null) html += '<br>';
      if (multiChain) html += '<span class="vi-seq-chain-lbl">Chain ' + r.chain + '</span><br>';
      prevChain = r.chain;
      chainPos  = 0;
    }}
    html += '<span class="vi-seq-aa" data-resi="' + r.resi + '" data-chain="' + r.chain
          + '" title="' + r.resn + ' ' + r.resi + '">' + r.aa + '</span>';
    chainPos++;
    if      (chainPos % 100 === 0) html += '<br>';
    else if (chainPos % 50  === 0) html += '<span class="vi-seq-sep">  </span>';
    else if (chainPos % 10  === 0) html += '<span class="vi-seq-sep"> </span>';
  }});
  document.getElementById('vi-seq-letters').innerHTML = html;
}})();

document.addEventListener('click', function(e) {{
  var el = e.target.closest('.vi-seq-aa');
  if (!el || !viewer) return;
  viSelectedResi  = parseInt(el.dataset.resi);
  viSelectedChain = el.dataset.chain;
  viColorSeq();
  viewer.center({{resi: viSelectedResi, chain: viSelectedChain}}, 500, false);
  viewer.render();
}});

window.addEventListener("DOMContentLoaded", function() {{
  viewer = $3Dmol.createViewer(document.getElementById("viewer"),
                               {{backgroundColor: "white", antialias: true}});
  viewer.addModel(PDB_DATA, "pdb");
  updateScene();
  viewer.zoomTo({{chain: LIG_CHAIN, resi: LIG_RESI}});
  viewer.render();
  document.querySelectorAll('input[name="vi-style"]').forEach(function(r) {{
    r.addEventListener('change', updateScene);
  }});
  document.querySelectorAll('input[name="vi-color"]').forEach(function(r) {{
    r.addEventListener('change', updateScene);
  }});
  document.getElementById('vi-spin').addEventListener('change', function() {{
    this.checked ? viewer.spin('y', 1) : viewer.spin(false);
  }});
  document.getElementById('vi-ortho').addEventListener('change', function() {{
    viewer.setProjection(this.checked ? 'orthographic' : 'perspective');
    viewer.render();
  }});
}});
</script>
</body>
</html>"""


# ─── 2D interaction fingerprint (matplotlib) ──────────────────────────────────

def plot_interaction_summary(report: LigandReport,
                              output_path: str | Path,
                              types: list[str] | None = None) -> Path:
    """
    Save a 2D interaction fingerprint chart for a LigandReport.

    Shows a dot-plot matrix (residues × interaction types) alongside a
    bar chart summarising total counts per type.

    Parameters
    ----------
    report :
        LigandReport from analyze_ligand().
    output_path :
        Output PNG (or any matplotlib-supported format) path.
    types :
        Subset of ["contact","hbond","pi","salt"] to include.  None = all.

    Returns
    -------
    Path to the saved image.
    """
    lig_key = (report.ligand_chain, report.ligand_resi)
    all_types = ["contact", "hbond", "pi", "salt"]
    show_types = types if types is not None else all_types

    # ── Aggregate per residue ─────────────────────────────────────────────────
    data: dict[tuple, dict] = defaultdict(
        lambda: {"resn": "?", "contact": [], "hbond": [], "pi": [], "salt": []}
    )

    for c in report.contacts:
        key = (c.chain2, c.resi2)
        if key != lig_key:
            data[key]["resn"] = c.resn2
            data[key]["contact"].append(c.distance)

    for h in report.hbonds:
        for chain, resi, resn in [
            (h.donor_chain,    h.donor_resi,    h.donor_resn),
            (h.acceptor_chain, h.acceptor_resi, h.acceptor_resn),
        ]:
            key = (chain, resi)
            if key != lig_key:
                data[key]["resn"] = resn
                data[key]["hbond"].append(h.distance)

    for p in report.pi_interactions:
        for chain, resi, resn in [(p.chain1, p.resi1, p.resn1),
                                   (p.chain2, p.resi2, p.resn2)]:
            key = (chain, resi)
            if key != lig_key:
                data[key]["resn"] = resn
                data[key]["pi"].append(p.center_distance)

    for s in report.salt_bridges:
        for chain, resi, resn in [(s.cation_chain, s.cation_resi, s.cation_resn),
                                   (s.anion_chain,  s.anion_resi,  s.anion_resn)]:
            key = (chain, resi)
            if key != lig_key:
                data[key]["resn"] = resn
                data[key]["salt"].append(s.distance)

    if not data:
        print("  No interactions to plot.")
        return Path(output_path)

    # Sort by residue number; filter to rows that have any data for show_types
    residue_keys = sorted(
        [k for k in data if any(data[k][t] for t in show_types)],
        key=lambda k: k[1],
    )
    if not residue_keys:
        print("  No interactions to plot for the requested types.")
        return Path(output_path)

    row_labels = [f"{data[k]['resn']}{k[1]}:{k[0]}" for k in residue_keys]
    _type_meta = {
        "contact": ("Contacts",     _COLORS["contact"]["hex"]),
        "hbond":   ("H-bonds",      _COLORS["hbond"]["hex"]),
        "pi":      ("Pi",           _COLORS["pi"]["hex"]),
        "salt":    ("Salt bridges", _COLORS["salt_cat"]["hex"]),
    }
    itypes     = show_types
    col_labels = [_type_meta[t][0] for t in itypes]
    dot_colors = [_type_meta[t][1] for t in itypes]

    n_res = len(residue_keys)
    fig_h = max(5, n_res * 0.38 + 2.5)
    fig, (ax_dot, ax_bar) = plt.subplots(
        1, 2, figsize=(11, fig_h),
        gridspec_kw={"width_ratios": [4, 1.2]},
    )

    # ── Dot plot ──────────────────────────────────────────────────────────────
    for yi, key in enumerate(residue_keys):
        row = data[key]
        for xi, itype in enumerate(itypes):
            dists = row[itype]
            if not dists:
                continue
            count = len(dists)
            min_d = min(dists)
            # Larger dot = more contacts; darker = closer
            size = 80 + count * 60
            alpha = min(1.0, max(0.5, 1.0 - (min_d - 2.0) / 4.0))
            ax_dot.scatter(xi, yi, s=size, c=dot_colors[xi],
                           alpha=alpha, zorder=4,
                           edgecolors="white", linewidths=0.8)
            if count > 1:
                ax_dot.text(xi, yi, str(count), ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold")
            # Show min distance below dot
            ax_dot.text(xi, yi - 0.33, f"{min_d:.1f}Å",
                        ha="center", va="center", fontsize=6.5,
                        color="#555555")

    ax_dot.set_yticks(range(n_res))
    ax_dot.set_yticklabels(row_labels, fontsize=9, fontfamily="monospace")
    ax_dot.set_xticks(range(len(itypes)))
    ax_dot.set_xticklabels(col_labels, fontsize=10, fontweight="bold")
    ax_dot.set_xlim(-0.6, len(itypes) - 0.4)
    ax_dot.set_ylim(-0.8, n_res - 0.2)
    ax_dot.grid(axis="x", color="#dddddd", linewidth=0.8)
    ax_dot.set_axisbelow(True)
    ax_dot.spines[["top", "right"]].set_visible(False)
    ax_dot.set_title(
        f"Interaction Fingerprint\n"
        f"{report.ligand_resn} {report.ligand_chain}:{report.ligand_resi}",
        fontsize=12, fontweight="bold", pad=10,
    )

    # Legend
    patches = [mpatches.Patch(color=dot_colors[i], label=col_labels[i])
               for i in range(len(itypes))]
    ax_dot.legend(handles=patches, loc="lower right", fontsize=8,
                  framealpha=0.9, edgecolor="#cccccc")

    # ── Summary bar chart ─────────────────────────────────────────────────────
    totals = {t: sum(len(data[k][t]) for k in residue_keys) for t in itypes}
    bars = ax_bar.barh(col_labels, [totals[t] for t in itypes],
                       color=dot_colors, edgecolor="white", linewidth=0.8,
                       height=0.55)
    ax_bar.bar_label(bars, fontsize=9, padding=4)
    max_val = max(totals.values()) if max(totals.values()) > 0 else 1
    ax_bar.set_xlim(0, max_val * 1.4)
    ax_bar.set_xlabel("Count", fontsize=9)
    ax_bar.set_title("Summary", fontsize=10, fontweight="bold")
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.tick_params(axis="y", labelsize=9)

    plt.tight_layout(pad=1.5)
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved fingerprint  → {output_path}")
    return output_path


# ─── LIGPLOT-style 2D diagram ─────────────────────────────────────────────────

def plot_ligand_2d(structure_path: str | Path,
                   report: LigandReport,
                   output_path: str | Path,
                   structure: Structure | None = None,
                   show_types: set[str] | None = None) -> Path:
    """
    Save a LIGPLOT-style 2D interaction diagram.

    Ligand heavy atoms are projected onto their best-fit plane via SVD and
    drawn as a 2D stick diagram with CPK atom colouring.  Interacting protein
    residues are placed evenly around the periphery as labelled colour-coded
    boxes.  Contacts use the classic LIGPLOT spoked-arc motif; H-bonds, pi
    interactions, and salt bridges are drawn as dashed lines.

    Parameters
    ----------
    structure_path :
        Path to the PDB or CIF file.
    report :
        LigandReport from analyze_ligand().
    output_path :
        Output PNG file path.
    structure :
        Pre-loaded BioPython Structure (optional).

    Returns
    -------
    Path to the saved PNG.
    """
    from matplotlib.lines import Line2D

    structure_path = Path(structure_path)
    if structure is None:
        structure = load_structure(structure_path)

    # ── Collect ligand heavy atoms ─────────────────────────────────────────────
    lig_atoms: list[tuple[str, str, np.ndarray]] = []   # (element, name, coord)
    for model in structure:
        for ch in model:
            if ch.id != report.ligand_chain:
                continue
            for res in ch:
                if res.id[1] != report.ligand_resi:
                    continue
                for atom in res:
                    elem = (atom.element or atom.get_name()[0]).strip().upper()
                    if elem and elem != "H":
                        lig_atoms.append((elem, atom.get_name().strip(),
                                          atom.get_vector().get_array()))
        break  # first model only

    if not lig_atoms:
        print("  No ligand heavy atoms found for 2D diagram.")
        return Path(output_path)

    # ── Project to 2D via SVD ──────────────────────────────────────────────────
    coords_3d   = np.array([a[2] for a in lig_atoms])
    centroid_3d = coords_3d.mean(axis=0)
    centered    = coords_3d - centroid_3d

    if centered.shape[0] >= 3:
        _, _, Vt = np.linalg.svd(centered)
        ax1, ax2 = Vt[0], Vt[1]
    else:
        ax1, ax2 = np.array([1, 0, 0]), np.array([0, 1, 0])

    lig_2d = np.column_stack([centered @ ax1, centered @ ax2])
    scale  = max(np.linalg.norm(lig_2d, axis=1).max(), 0.1)
    lig_2d /= scale          # normalise to ~unit circle

    # ── Infer covalent bonds (<1.9 Å between heavy atoms in 3D) ──────────────
    bonds: list[tuple[int, int]] = []
    for i in range(len(lig_atoms)):
        for j in range(i + 1, len(lig_atoms)):
            if np.linalg.norm(coords_3d[i] - coords_3d[j]) < 1.9:
                bonds.append((i, j))

    # ── Atom name -> 2D position lookup ──────────────────────────────────────
    name_to_2d: dict[str, np.ndarray] = {
        name: lig_2d[i] for i, (_, name, _) in enumerate(lig_atoms)
    }

    # ── Collect interacting residues, types, and specific ligand atom ─────────
    lig_key = (report.ligand_chain, report.ligand_resi)
    _IPRI = {"hbond": 0, "salt": 1, "pi": 2, "contact": 3}

    res_types:      dict[tuple, set[str]] = defaultdict(set)
    res_names:      dict[tuple, str]      = {}
    res_lig_atom:   dict[tuple, str]      = {}   # best-priority ligand atom per residue
    res_best_itype: dict[tuple, str]      = {}

    def _add(key, resn, itype, lig_atom):
        if show_types is not None and itype not in show_types:
            return
        res_types[key].add(itype)
        res_names[key] = resn
        if _IPRI[itype] < _IPRI.get(res_best_itype.get(key, "contact"), 3):
            res_best_itype[key] = itype
            res_lig_atom[key]   = lig_atom.strip()

    for c in report.contacts:
        if (c.chain1, c.resi1) == lig_key:
            _add((c.chain2, c.resi2), c.resn2, "contact", c.atom1)
        elif (c.chain2, c.resi2) == lig_key:
            _add((c.chain1, c.resi1), c.resn1, "contact", c.atom2)

    for h in report.hbonds:
        if (h.donor_chain, h.donor_resi) == lig_key:
            _add((h.acceptor_chain, h.acceptor_resi), h.acceptor_resn, "hbond", h.donor_atom)
        elif (h.acceptor_chain, h.acceptor_resi) == lig_key:
            _add((h.donor_chain, h.donor_resi), h.donor_resn, "hbond", h.acceptor_atom)

    for p in report.pi_interactions:
        if (p.chain1, p.resi1) == lig_key:
            _add((p.chain2, p.resi2), p.resn2, "pi", p.ring1_label.split(",")[0])
        elif (p.chain2, p.resi2) == lig_key:
            _add((p.chain1, p.resi1), p.resn1, "pi", p.ring2_label.split(",")[0])

    for s in report.salt_bridges:
        if (s.cation_chain, s.cation_resi) == lig_key:
            _add((s.anion_chain, s.anion_resi), s.anion_resn, "salt", s.cation_atom)
        elif (s.anion_chain, s.anion_resi) == lig_key:
            _add((s.cation_chain, s.cation_resi), s.cation_resn, "salt", s.anion_atom)

    if not res_types:
        print("  No interactions found for 2D diagram.")
        return Path(output_path)

    # ── Place residues evenly around a circle ─────────────────────────────────
    res_keys = sorted(res_types.keys(), key=lambda k: k[1])
    n_res    = len(res_keys)
    RING_R   = 2.3
    angles   = np.linspace(0, 2 * np.pi, n_res, endpoint=False) - np.pi / 2
    res_pos: dict[tuple, np.ndarray] = {
        key: np.array([RING_R * np.cos(a), RING_R * np.sin(a)])
        for key, a in zip(res_keys, angles)
    }

    # ── Figure setup ──────────────────────────────────────────────────────────
    fig_size = min(max(9, 2.5 + n_res * 0.5), 18)
    fig, ax  = plt.subplots(figsize=(fig_size, fig_size))
    ax.set_aspect("equal");  ax.axis("off")
    lim = RING_R + 1.2
    ax.set_xlim(-lim, lim);  ax.set_ylim(-lim, lim)
    ax.set_title(
        f"LIGPLOT-style 2D Diagram\n"
        f"{report.ligand_resn} {report.ligand_chain}:{report.ligand_resi}",
        fontsize=12, fontweight="bold", pad=12,
    )

    # ── Draw interaction lines (behind everything) ────────────────────────────
    # Line colours match the 3D viewer palette
    _LINE_STYLE = {
        "hbond":   (_DASH_COLORS["hbond"], "--", 1.8),   # yellow
        "salt":    (_DASH_COLORS["salt"],  "--", 1.8),   # orange
        "pi":      (_DASH_COLORS["pi"],    "--", 1.6),   # magenta
        "contact": (_COLORS["contact"]["hex"], "-", 1.2), # cyan
    }
    for key in res_keys:
        types = res_types[key]
        rpos  = res_pos[key]

        # Pick highest-priority type for this residue
        itype = min(types, key=lambda t: _IPRI[t])
        color, ls, lw = _LINE_STYLE[itype]

        # Anchor at the specific interacting ligand atom (fallback: origin)
        atom_nm  = res_lig_atom.get(key, "")
        atom_pos = name_to_2d.get(atom_nm, np.array([0.0, 0.0]))

        direc  = atom_pos - rpos
        d_norm = np.linalg.norm(direc)
        if d_norm > 0:
            direc /= d_norm

        ax.plot([atom_pos[0], rpos[0]], [atom_pos[1], rpos[1]],
                ls, color=color, linewidth=lw, zorder=1, alpha=0.9)

        if types == {"contact"}:          # pure hydrophobic → spoked arc
            _draw_spoked_arc(ax, rpos, direc)

    # ── Draw ligand bonds (sticks) ────────────────────────────────────────────
    for i, j in bonds:
        ei = lig_atoms[i][0]; ej = lig_atoms[j][0]
        ci = _CPK_COLORS.get(ei, _CPK_COLORS.get(ei[:1], "#FF1493"))
        cj = _CPK_COLORS.get(ej, _CPK_COLORS.get(ej[:1], "#FF1493"))
        xi, yi = lig_2d[i]; xj, yj = lig_2d[j]
        mx, my = (xi+xj)/2, (yi+yj)/2
        # Draw each half in its atom's CPK colour
        ax.plot([xi, mx], [yi, my], color=ci, linewidth=6.0,
                solid_capstyle="round", zorder=2)
        ax.plot([mx, xj], [my, yj], color=cj, linewidth=6.0,
                solid_capstyle="round", zorder=2)

    # ── Heteroatom element labels (no circles) ────────────────────────────────
    for idx, (elem, _name, _) in enumerate(lig_atoms):
        if elem == "C":
            continue
        color = _CPK_COLORS.get(elem, _CPK_COLORS.get(elem[:1], "#FF1493"))
        x, y  = lig_2d[idx]
        ax.text(x, y, elem, ha="center", va="center",
                fontsize=7.0, color=color, zorder=4, fontweight="bold")

    # ── Draw residue boxes ────────────────────────────────────────────────────
    for key in res_keys:
        types = res_types[key]
        rpos  = res_pos[key]
        label = f"{res_names.get(key,'?')}{key[1]}:{key[0]}"
        if "hbond" in types:
            bg = _COLORS["hbond"]["hex"]
        elif "salt" in types:
            bg = _COLORS["salt_cat"]["hex"]
        elif "pi" in types:
            bg = _COLORS["pi"]["hex"]
        else:
            bg = "#777777"
        ax.text(rpos[0], rpos[1], label,
                ha="center", va="center", fontsize=8.5,
                fontfamily="monospace", fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=0.38", facecolor=bg,
                          edgecolor="white", linewidth=1.5),
                zorder=5)

    # ── Legend ─────────────────────────────────────────────────────────────────
    legend_items = [
        Line2D([0], [0], ls="--", color=_DASH_COLORS["hbond"],        lw=1.8, label="H-bond"),
        Line2D([0], [0], ls="--", color=_DASH_COLORS["pi"],           lw=1.6, label="Pi interaction"),
        Line2D([0], [0], ls="--", color=_DASH_COLORS["salt"],         lw=1.8, label="Salt bridge"),
        Line2D([0], [0], ls="-",  color=_COLORS["contact"]["hex"],    lw=1.2, label="Hydrophobic contact"),
    ]
    ax.legend(handles=legend_items, loc="lower right", fontsize=8.5,
              framealpha=0.9, edgecolor="#cccccc")

    plt.tight_layout()
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved 2D diagram   → {output_path}")
    return output_path


# ─── Per-interaction-type separate images ─────────────────────────────────────

def save_type_images(structure_path: str | Path,
                     report: LigandReport,
                     outdir: str | Path,
                     stem: str,
                     structure: Structure | None = None,
                     width: int = 900,
                     height: int = 700) -> list[dict]:
    """
    Generate per-interaction-type HTML viewers and filtered fingerprint PNGs.

    For each interaction type that has at least one detected interaction, writes:
      - ``{stem}_{type}.html``        — 3D viewer with only that type pre-checked
      - ``{stem}_{type}_fingerprint.png`` — fingerprint filtered to that type

    Parameters
    ----------
    structure_path :
        PDB or CIF file.
    report :
        LigandReport from analyze_ligand().
    outdir :
        Output directory (must already exist).
    stem :
        Filename stem (e.g. ``"ATP_A501"``).
    structure :
        Pre-loaded BioPython Structure (optional).
    width, height :
        HTML viewer dimensions.

    Returns
    -------
    List of dicts with keys ``"type"``, ``"html"``, ``"png"`` per type generated.
    """
    structure_path = Path(structure_path)
    outdir         = Path(outdir)
    if structure is None:
        structure = load_structure(structure_path)

    # Map HTML type-id → (interaction list, filter key for fingerprint)
    type_map = [
        ("contacts", report.contacts,        ["contact"]),
        ("hbonds",   report.hbonds,          ["hbond"]),
        ("pi",       report.pi_interactions, ["pi"]),
        ("salt",     report.salt_bridges,    ["salt"]),
    ]

    results = []
    for type_id, items, filter_keys in type_map:
        if not items:
            continue

        # Type-specific 3D HTML (only this checkbox pre-checked)
        html_path = outdir / f"{stem}_{type_id}.html"
        html = _render_html(structure_path, report, structure, width, height,
                            active_types={type_id})
        html_path.write_text(html)
        print(f"  Saved {type_id:<10} HTML → {html_path}")

        # Type-specific fingerprint PNG (only matching columns shown)
        png_path = outdir / f"{stem}_{type_id}_fingerprint.png"
        plot_interaction_summary(report, png_path, types=filter_keys)

        results.append({"type": type_id, "html": html_path, "png": png_path})

    return results


# ─── High-level entry point ────────────────────────────────────────────────────

def visualize(structure_path: str | Path,
              ligand_sel: str | None = None,
              protein_sel: str = "polymer",
              outdir: str | Path = ".",
              prefix: str = "",
              width: int = 900,
              height: int = 700,
              contact_cutoff: float = 4.5) -> list[dict]:
    """
    Run full interaction analysis and generate both HTML and PNG outputs.

    Parameters
    ----------
    structure_path :
        PDB or CIF file.
    ligand_sel :
        PyMOL-like selection for the ligand(s).  If ``None``, auto-detects
        primary ligands (non-excipient, ≥7 heavy atoms) sorted largest first.
    protein_sel :
        PyMOL-like selection for the binding partner. Default ``"polymer"``.
    outdir :
        Output directory (created if missing).
    prefix :
        Optional filename prefix.
    width, height :
        HTML viewer dimensions.
    contact_cutoff :
        Contact distance cutoff in Å.

    Returns
    -------
    List of dicts with keys ``"report"``, ``"html"``, ``"png"``, ``"diagram2d"``,
    ``"type_images"`` per ligand.
    """
    structure_path = Path(structure_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    structure = load_structure(structure_path)

    # ── Auto-detect primary ligands when no selection is given ────────────────
    if ligand_sel is None:
        candidates = primary_ligands(structure)
        if not candidates:
            print("No primary ligands found (all HETATM residues are water, "
                  "excipients, or < 7 heavy atoms).  Use --analyze to specify "
                  "a ligand explicitly.")
            return []
        # Build a precise selection covering all detected primary ligands
        parts = [f"(resn {r.get_resname().strip()} and chain {r.get_parent().id} "
                 f"and resi {r.id[1]})" for r in candidates]
        ligand_sel = " or ".join(parts)
        names = ", ".join(
            f"{r.get_resname().strip()} {r.get_parent().id}:{r.id[1]}"
            for r in candidates
        )
        print(f"Auto-detected primary ligand(s): {names}")

    reports = analyze_ligand(structure, ligand_sel, protein_sel=protein_sel,
                             contact_cutoff=contact_cutoff,
                             pdb_path=structure_path)

    if not reports:
        print(f"No residues matched: {ligand_sel!r}")
        return []

    results = []
    for report in reports:
        tag  = f"{report.ligand_resn}_{report.ligand_chain}{report.ligand_resi}"
        stem = f"{prefix}{tag}" if prefix else tag

        print(f"\n{report.summary()}")

        html_path = save_html(structure_path, report, outdir / f"{stem}.html",
                              structure=structure, width=width, height=height)
        png_path  = plot_interaction_summary(report, outdir / f"{stem}_fingerprint.png")
        diag2d    = plot_ligand_2d(structure_path, report,
                                   outdir / f"{stem}_2d.png", structure=structure)
        type_imgs = save_type_images(structure_path, report, outdir, stem,
                                     structure=structure, width=width, height=height)

        results.append({"report": report, "html": html_path, "png": png_path,
                        "diagram2d": diag2d, "type_images": type_imgs})

    return results


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate 3D HTML and 2D fingerprint images for structural interactions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 visualize_interactions.py 1abc.pdb --analyze "organic"
  python3 visualize_interactions.py 1abc.pdb --analyze "resn ATP" --outdir ./images
  python3 visualize_interactions.py 1abc.pdb --analyze "resn ATP and chain B" --width 1200

PyMOL selection syntax is supported (see analyze_ligands.py for full reference).
""",
    )

    parser.add_argument("structure",  help="PDB or mmCIF file")
    parser.add_argument("--analyze",  default=None, metavar="SELECTION",
                        help="Ligand selection, e.g. \"resn ATP\" or \"organic\". "
                             "If omitted, auto-detects primary ligands (non-excipient, "
                             "≥7 heavy atoms) sorted largest first.")
    parser.add_argument("--protein",  default="polymer",
                        help="Binding partner selection (default: polymer)")
    parser.add_argument("--cutoff",   type=float, default=4.5,
                        help="Contact distance cutoff in Å (default: 4.5)")
    parser.add_argument("--outdir",   default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("--prefix",   default="",
                        help="Optional filename prefix")
    parser.add_argument("--width",    type=int, default=900,
                        help="HTML viewer width in pixels (default: 900)")
    parser.add_argument("--height",   type=int, default=700,
                        help="HTML viewer height in pixels (default: 700)")

    args = parser.parse_args()

    results = visualize(
        structure_path  = args.structure,
        ligand_sel      = args.analyze,   # None triggers auto-detection
        protein_sel     = args.protein,
        outdir          = args.outdir,
        prefix          = args.prefix,
        width           = args.width,
        height          = args.height,
        contact_cutoff  = args.cutoff,
    )

    if results:
        n_type = sum(len(r["type_images"]) for r in results)
        n_files = len(results) * 3 + n_type * 2   # html + fingerprint + 2d + per-type
        print(f"\nGenerated {n_files} files in '{args.outdir}'")


if __name__ == "__main__":
    main()
