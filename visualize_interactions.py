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
    SelectionParser, analyze_ligand, load_structure,
    WATER_NAMES,
)
from Bio.PDB.Structure import Structure

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
              height: int = 700) -> Path:
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
    html = _render_html(structure_path, report, structure, width, height)
    output_path = Path(output_path)
    output_path.write_text(html)
    print(f"  Saved 3D viewer    → {output_path}")
    return output_path


def _render_html(structure_path: Path, report: LigandReport,
                 structure: Structure, width: int, height: int) -> str:
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

    # ── Header subtitle ───────────────────────────────────────────────────────
    subtitle = (f"{report.ligand_resn} {report.ligand_chain}:{report.ligand_resi}"
                f" &nbsp;·&nbsp; {report.n_contacts} contacts"
                f" &nbsp;·&nbsp; {report.n_hbonds} H-bonds"
                f" &nbsp;·&nbsp; {report.n_pi} pi"
                f" &nbsp;·&nbsp; {report.n_salt_bridges} salt bridges")

    # ── Serialise data for JS ─────────────────────────────────────────────────
    J = json.dumps
    water_js = J(sorted(WATER_NAMES))

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
</style>
</head>
<body>
<div class="header">
  <h2>Interaction Viewer</h2>
  <div class="subtitle">{subtitle}</div>
</div>
<div class="body">
  <div id="viewer-wrap"><div id="viewer"></div></div>
  <div class="panel">
    <div class="toggles">
      <div class="toggles-title">Show / Hide</div>
      <label class="tog-row">
        <input type="checkbox" id="tog-contacts" checked onchange="updateScene()">
        <span class="tog-label">
          <span class="tog-dot" style="background:{_COLORS['contact']['hex']}"></span>
          Contacts ({len(contact_items)})
        </span>
      </label>
      <label class="tog-row">
        <input type="checkbox" id="tog-hbonds" checked onchange="updateScene()">
        <span class="tog-label">
          <span class="tog-dot" style="background:{_COLORS['hbond']['hex']}"></span>
          H-bonds ({len(hbond_items)})
        </span>
      </label>
      <label class="tog-row">
        <input type="checkbox" id="tog-pi" checked onchange="updateScene()">
        <span class="tog-label">
          <span class="tog-dot" style="background:{_COLORS['pi']['hex']}"></span>
          Pi interactions ({len(pi_items)})
        </span>
      </label>
      <label class="tog-row">
        <input type="checkbox" id="tog-salt" checked onchange="updateScene()">
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

var viewer;

function updateScene() {{
  // 1. Reset: base cartoon
  viewer.setStyle({{}}, {{cartoon: {{color: "#b8b8b8", opacity: 0.75}}}});
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
}}

function setAll(on) {{
  ["contacts","hbonds","pi","salt"].forEach(function(t) {{
    document.getElementById("tog-"+t).checked = on;
  }});
  updateScene();
}}

window.addEventListener("DOMContentLoaded", function() {{
  viewer = $3Dmol.createViewer(document.getElementById("viewer"),
                               {{backgroundColor: "white"}});
  viewer.addModel(PDB_DATA, "pdb");
  updateScene();
  viewer.zoomTo({{chain: LIG_CHAIN, resi: LIG_RESI}});
  viewer.render();
}});
</script>
</body>
</html>"""


# ─── 2D interaction fingerprint (matplotlib) ──────────────────────────────────

def plot_interaction_summary(report: LigandReport,
                              output_path: str | Path) -> Path:
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

    Returns
    -------
    Path to the saved image.
    """
    lig_key = (report.ligand_chain, report.ligand_resi)

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

    # Sort by residue number
    residue_keys = sorted(data.keys(), key=lambda k: k[1])
    row_labels = [f"{data[k]['resn']}{k[1]}:{k[0]}" for k in residue_keys]
    itypes = ["contact", "hbond", "pi", "salt"]
    col_labels = ["Contacts", "H-bonds", "Pi", "Salt bridges"]
    dot_colors = [_COLORS["contact"]["hex"], _COLORS["hbond"]["hex"],
                  _COLORS["pi"]["hex"],      _COLORS["salt_cat"]["hex"]]

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
            alpha = max(0.5, 1.0 - (min_d - 2.0) / 4.0)
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


# ─── High-level entry point ────────────────────────────────────────────────────

def visualize(structure_path: str | Path,
              ligand_sel: str,
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
        PyMOL-like selection for the ligand(s).
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
    List of dicts with keys ``"report"``, ``"html"``, ``"png"`` per ligand.
    """
    structure_path = Path(structure_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    structure = load_structure(structure_path)
    reports = analyze_ligand(structure, ligand_sel, protein_sel=protein_sel,
                             contact_cutoff=contact_cutoff,
                             pdb_path=structure_path)

    if not reports:
        print(f"No residues matched: {ligand_sel!r}")
        return []

    results = []
    for report in reports:
        tag = f"{report.ligand_resn}_{report.ligand_chain}{report.ligand_resi}"
        stem = f"{prefix}{tag}" if prefix else tag

        print(f"\n{report.summary()}")

        html_path = save_html(structure_path, report, outdir / f"{stem}.html",
                              structure=structure, width=width, height=height)
        png_path  = plot_interaction_summary(report, outdir / f"{stem}_fingerprint.png")

        results.append({"report": report, "html": html_path, "png": png_path})

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
    parser.add_argument("--analyze",  required=True, metavar="SELECTION",
                        help="Ligand selection, e.g. \"resn ATP\" or \"organic\"")
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
        ligand_sel      = args.analyze,
        protein_sel     = args.protein,
        outdir          = args.outdir,
        prefix          = args.prefix,
        width           = args.width,
        height          = args.height,
        contact_cutoff  = args.cutoff,
    )

    if results:
        print(f"\nGenerated {len(results) * 2} files in '{args.outdir}'")


if __name__ == "__main__":
    main()
