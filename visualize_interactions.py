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
import py3Dmol

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
    """Simulate a dashed line with alternating cylinders."""
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


# ─── 3D HTML viewer (py3Dmol) ─────────────────────────────────────────────────

def build_view(structure_path: str | Path,
               report: LigandReport,
               structure: Structure | None = None,
               width: int = 900,
               height: int = 700) -> py3Dmol.view:
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
    view = build_view(structure_path, report, structure, width, height)
    html = _render_html(view, report, width, height)
    output_path = Path(output_path)
    output_path.write_text(html)
    print(f"  Saved 3D viewer    → {output_path}")
    return output_path


def _render_html(view: py3Dmol.view, report: LigandReport,
                 width: int, height: int) -> str:
    """Wrap py3Dmol view in a full HTML page with a legend."""
    legend_items = [
        ("Ligand",        _COLORS["ligand"]["hex"]),
        ("Contacts",      _COLORS["contact"]["hex"]),
        ("H-bonds",       _COLORS["hbond"]["hex"]),
        ("Pi interaction",_COLORS["pi"]["hex"]),
        ("Salt bridge (+)",_COLORS["salt_cat"]["hex"]),
        ("Salt bridge (−)",_COLORS["salt_ani"]["hex"]),
    ]
    legend_html = "".join(
        f'<span style="display:inline-block;width:14px;height:14px;'
        f'background:{c};border-radius:3px;margin-right:5px;vertical-align:middle"></span>'
        f'<span style="margin-right:16px;font-size:13px">{label}</span>'
        for label, c in legend_items
    )

    title = (f"{report.ligand_resn} {report.ligand_chain}:{report.ligand_resi} — "
             f"{report.n_contacts} contacts · {report.n_hbonds} H-bonds · "
             f"{report.n_pi} pi · {report.n_salt_bridges} salt bridges")

    # Extract the inner viewer HTML from py3Dmol
    inner = view._make_html()

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{report.ligand_resn} Interactions</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 0; background: #f5f5f5; }}
  .header {{ background: #2c3e50; color: white; padding: 12px 20px; }}
  .header h2 {{ margin: 0 0 6px 0; font-size: 16px; }}
  .legend {{ padding: 8px 20px; background: white; border-bottom: 1px solid #ddd; }}
  .viewer {{ display: flex; justify-content: center; padding: 10px; }}
  iframe {{ border: none; border-radius: 4px; box-shadow: 0 2px 8px rgba(0,0,0,0.15); }}
</style>
</head>
<body>
<div class="header">
  <h2>Interaction Viewer</h2>
  <div style="font-size:13px;opacity:0.9">{title}</div>
</div>
<div class="legend">{legend_html}</div>
<div class="viewer">
  <iframe srcdoc="{_escape(inner)}" width="{width}" height="{height}"></iframe>
</div>
</body>
</html>"""


def _escape(s: str) -> str:
    return s.replace("&", "&amp;").replace('"', "&quot;")


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
