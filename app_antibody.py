#!/usr/bin/env python3
"""
app_antibody.py — Antibody Structure Analyzer.

Extends the PDB Structure Analyzer with antibody-specific analyses:
  • Heuristic VH/VL chain classification (no HMMER required)
  • CDR identification using conserved anchor residues (Chothia approximate)
  • CDR-colored 3D viewer
  • VH–VL interface contacts
  • Paratope–epitope contacts (when an antigen chain is present)

CDR boundaries are estimated from conserved C and W anchor positions plus
canonical offsets; results are approximate (±2 residues on CDR edges).

Run:
    streamlit run app_antibody.py
"""

from __future__ import annotations

import json
import re
import sys
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from Bio.PDB import PDBParser, NeighborSearch

TOOLKIT_DIR = Path(__file__).parent
sys.path.insert(0, str(TOOLKIT_DIR))

from summarize_structures import (
    bfactor_stats,
    buried_surface_areas,
    extract_fasta,
    find_missing_residues,
    load_structure as ss_load,
    plot_bfactor,
    plot_bsa_matrix,
    plot_ramachandran,
    radius_of_gyration,
    ramachandran_data,
    summarize_structure,
)
from analyze_ligands import (
    analyze_ligand,
    get_ligands,
    load_structure as al_load,
    primary_ligands,
)
from visualize_interactions import (
    _render_html,
    plot_interaction_summary,
    plot_ligand_2d,
)

# ── Amino acid mapping ─────────────────────────────────────────────────────────
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "MSE": "M",
}

# ── Color palettes ─────────────────────────────────────────────────────────────
_CHAIN_COLORS = {
    "VH": "#4e79a7",
    "VL": "#f28e2b",
    "Antigen": "#e15759",
    "Unknown": "#888888",
}
_CDR_COLORS = {
    "CDR-H1": "#ff4444",   # red
    "CDR-H2": "#ff8800",   # orange
    "CDR-H3": "#ffee00",   # yellow
    "CDR-L1": "#44cc44",   # green
    "CDR-L2": "#4488ff",   # blue
    "CDR-L3": "#cc44ff",   # violet
}

st.markdown("""
<style>
.stTabs [data-baseweb="tab-list"] { gap: 6px; }
.stTabs [data-baseweb="tab"] { padding: 8px 18px; border-radius: 6px 6px 0 0; }
</style>
""", unsafe_allow_html=True)

_NO_PDB = "Upload a PDB file above to get started."


# ═══════════════════════════════════════════════════════════════════════════════
# Heuristic CDR detection
# ═══════════════════════════════════════════════════════════════════════════════

def _classify_antibody_chain(seq: str) -> dict:
    """
    Heuristic VH/VL classification and CDR boundary detection.

    Anchors used:
      - FR2 tryptophan motif  → chain type (VH: W-[VIL]-[RK], VL: W-[YF]-QQ)
      - C1 ≈ Kabat C22/23    → ~14 residues before FR2-W
      - C2 ≈ Kabat C92/88    → ~69 (VH) / 65 (VL) residues after C1

    CDR offsets from anchors follow Chothia/Kabat canonical ranges.
    Returns dict with keys: type, c1_pos, c2_pos, cdr1, cdr2, cdr3
    where cdrN is a (start, end) half-open slice, or None.
    """
    s = seq.upper()
    n = len(s)
    result = {"type": "Unknown", "c1_pos": None, "c2_pos": None,
              "cdr1": None, "cdr2": None, "cdr3": None}

    if n < 70:
        return result

    # Classify by FR2 tryptophan context
    vh_fr2 = re.search(r'W[VILMF][RK][QKE]', s)
    vl_fr2 = re.search(r'W[YF]QQ', s)

    if vh_fr2 and vl_fr2:
        chain_type = "VH" if vh_fr2.start() <= vl_fr2.start() else "VL"
        fr2_m = vh_fr2 if chain_type == "VH" else vl_fr2
    elif vh_fr2:
        chain_type, fr2_m = "VH", vh_fr2
    elif vl_fr2:
        chain_type, fr2_m = "VL", vl_fr2
    else:
        return result

    result["type"] = chain_type
    w36 = fr2_m.start()

    # C1 (≈ Kabat C22/23): ~14 (VH) or ~13 (VL) residues before FR2-W
    c1_off = 14 if chain_type == "VH" else 13
    lo = max(0, w36 - c1_off - 4)
    hi = min(n, w36 - c1_off + 4)
    c1_rel = s[lo:hi].rfind("C")
    c1 = (lo + c1_rel) if c1_rel >= 0 else max(0, w36 - c1_off)
    result["c1_pos"] = c1

    # C2 (≈ Kabat C92/88): ~69 (VH) or ~65 (VL) residues after C1
    c2_off = 69 if chain_type == "VH" else 65
    lo2 = max(0, c1 + c2_off - 5)
    hi2 = min(n, c1 + c2_off + 5)
    c2_rel = s[lo2:hi2].find("C")
    c2 = (lo2 + c2_rel) if c2_rel >= 0 else min(n - 1, c1 + c2_off)
    result["c2_pos"] = c2

    def _clamp(s_: int, e_: int) -> tuple[int, int] | None:
        s_, e_ = max(0, min(s_, n)), max(0, min(e_, n))
        return (s_, e_) if s_ < e_ else None

    if chain_type == "VH":
        result["cdr1"] = _clamp(c1 + 4, min(w36 - 3, c1 + 13))     # H1 Chothia 26-32
        result["cdr2"] = _clamp(w36 + 16, min(c2 - 12, w36 + 26))   # H2 Kabat 52-58
        m = re.search(r'WG[QA][GE]T', s[c2:min(c2 + 25, n)])
        cdr3_e = (c2 + m.start()) if m else min(c2 + 16, n)
        result["cdr3"] = _clamp(c2 + 3, cdr3_e)                      # H3 Kabat 95-102
    else:
        result["cdr1"] = _clamp(c1 + 4, min(w36 - 3, c1 + 18))      # L1 Kabat 24-34
        result["cdr2"] = _clamp(w36 + 14, min(c2 - 12, w36 + 21))   # L2 Kabat 50-56
        m = re.search(r'FG[QA][GE]T', s[c2:min(c2 + 20, n)])
        cdr3_e = (c2 + m.start()) if m else min(c2 + 12, n)
        result["cdr3"] = _clamp(c2 + 3, cdr3_e)                      # L3 Kabat 89-97

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# Cached analysis functions (structure + ligands — same as app_structure.py)
# ═══════════════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def _run_structure_analysis(data: bytes) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data); p = Path(f.name)
    try:
        ss_structure = ss_load(str(p))
        return {
            "summary":      summarize_structure(p),
            "bfactor":      bfactor_stats(ss_structure),
            "ramachandran": ramachandran_data(ss_structure),
            "rg":           radius_of_gyration(ss_structure),
            "fasta":        extract_fasta(p),
            "missing":      find_missing_residues(p),
        }
    finally:
        p.unlink(missing_ok=True)


@st.cache_data(show_spinner=False)
def _run_ligand_analysis(data: bytes, cutoff: float):
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data); p = Path(f.name)
    try:
        structure = al_load(str(p))
        candidates = primary_ligands(structure)
        if not candidates:
            return []
        parts = [
            f"(resn {r.get_resname().strip()} and chain {r.get_parent().id} "
            f"and resi {r.id[1]})"
            for r in candidates
        ]
        return analyze_ligand(structure, " or ".join(parts), pdb_path=p, contact_cutoff=cutoff)
    finally:
        p.unlink(missing_ok=True)


@st.cache_data(show_spinner=False)
def _run_bsa(data: bytes, probe_radius: float, n_points: int):
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data); p = Path(f.name)
    try:
        structure = ss_load(str(p))
        chain_ids = [ch.id for ch in structure[0]]
        results   = buried_surface_areas(structure, probe_radius=probe_radius, n_points=n_points)
        return results, chain_ids
    finally:
        p.unlink(missing_ok=True)


@st.cache_data(show_spinner=False)
def _list_all_ligands(data: bytes) -> list[dict]:
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data); p = Path(f.name)
    try:
        structure = al_load(str(p))
        ligs = get_ligands(structure, exclude_water=True)
        return [{"resn": res.get_resname().strip(), "chain": res.get_parent().id,
                 "resi": res.id[1], "atoms": sum(1 for _ in res.get_atoms())}
                for res in ligs]
    finally:
        p.unlink(missing_ok=True)


# ── Antibody-specific analysis ─────────────────────────────────────────────────

@st.cache_data(show_spinner=False)
def _run_antibody_analysis(data: bytes) -> dict:
    """Classify chains, detect CDRs, find VH-VL and paratope-epitope contacts."""
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data); p = Path(f.name)
    try:
        parser = PDBParser(QUIET=True)
        struct  = parser.get_structure("ab", str(p))
        model   = struct[0]

        chain_info: dict[str, dict] = {}
        for chain in model:
            residues = [
                (res.get_id()[1], res.get_resname(), _AA3TO1[res.get_resname()])
                for res in chain
                if res.get_id()[0] == " " and res.get_resname() in _AA3TO1
            ]
            if len(residues) < 30:
                continue
            seq = "".join(r[2] for r in residues)
            ab  = _classify_antibody_chain(seq)

            entry: dict = {
                "chain_id": chain.id,
                "length":   len(residues),
                "sequence": seq,
                "residues": residues,
                "type":     ab["type"],
            }
            for cdr_key in ("cdr1", "cdr2", "cdr3"):
                rng = ab[cdr_key]
                if rng:
                    s_, e_ = rng
                    entry[f"{cdr_key}_seq"]   = seq[s_:e_]
                    entry[f"{cdr_key}_resis"]  = [residues[i][0] for i in range(s_, e_) if i < len(residues)]
                else:
                    entry[f"{cdr_key}_seq"]   = None
                    entry[f"{cdr_key}_resis"]  = []
            chain_info[chain.id] = entry

        vh_chains = [cid for cid, v in chain_info.items() if v["type"] == "VH"]
        vl_chains = [cid for cid, v in chain_info.items() if v["type"] == "VL"]
        ag_chains = [cid for cid, v in chain_info.items() if v["type"] == "Unknown"]

        # ── VH-VL interface contacts ───────────────────────────────────────────
        vhvl_contacts: list[dict] = []
        if vh_chains and vl_chains:
            vh_id, vl_id = vh_chains[0], vl_chains[0]
            vh_atoms = list(model[vh_id].get_atoms())
            vl_atoms = list(model[vl_id].get_atoms())
            ns   = NeighborSearch(vh_atoms + vl_atoms)
            seen: set = set()
            for atom in vh_atoms:
                for n_atom in ns.search(atom.coord, 4.5):
                    n_res = n_atom.get_parent()
                    if n_res.get_parent().id == vl_id:
                        a_res = atom.get_parent()
                        key   = (a_res.get_id()[1], n_res.get_id()[1])
                        if key not in seen:
                            seen.add(key)
                            vhvl_contacts.append({
                                "VH chain": vh_id, "VH resi": a_res.get_id()[1], "VH resn": a_res.get_resname(),
                                "VL chain": vl_id, "VL resi": n_res.get_id()[1], "VL resn": n_res.get_resname(),
                                "dist (Å)": round(float(atom - n_atom), 2),
                            })
            vhvl_contacts.sort(key=lambda x: x["dist (Å)"])

        # ── Paratope-epitope contacts ──────────────────────────────────────────
        paratope_contacts: list[dict] = []
        ab_ids = vh_chains + vl_chains
        if ab_ids and ag_chains:
            ab_atoms = [a for cid in ab_ids for a in model[cid].get_atoms()]
            ag_atoms = [a for cid in ag_chains for a in model[cid].get_atoms()]
            cdr_resi_map = {
                cid: set(
                    resi
                    for k in ("cdr1_resis", "cdr2_resis", "cdr3_resis")
                    for resi in chain_info[cid].get(k, [])
                )
                for cid in ab_ids
            }
            ns   = NeighborSearch(ab_atoms + ag_atoms)
            seen = set()
            for atom in ag_atoms:
                for n_atom in ns.search(atom.coord, 4.5):
                    n_res      = n_atom.get_parent()
                    n_chain_id = n_res.get_parent().id
                    if n_chain_id in ab_ids:
                        ag_res   = atom.get_parent()
                        ag_cid   = ag_res.get_parent().id
                        key      = (ag_cid, ag_res.get_id()[1], n_chain_id, n_res.get_id()[1])
                        if key not in seen:
                            seen.add(key)
                            in_cdr = n_res.get_id()[1] in cdr_resi_map.get(n_chain_id, set())
                            paratope_contacts.append({
                                "Ag chain": ag_cid, "Ag resi": ag_res.get_id()[1], "Ag resn": ag_res.get_resname(),
                                "Ab chain": n_chain_id, "Ab type": chain_info[n_chain_id]["type"],
                                "Ab resi": n_res.get_id()[1], "Ab resn": n_res.get_resname(),
                                "In CDR": "✓" if in_cdr else "",
                                "dist (Å)": round(float(atom - n_atom), 2),
                            })
            paratope_contacts.sort(key=lambda x: x["dist (Å)"])

        return {
            "chains":            chain_info,
            "vh_chains":         vh_chains,
            "vl_chains":         vl_chains,
            "ag_chains":         ag_chains,
            "vhvl_contacts":     vhvl_contacts,
            "paratope_contacts": paratope_contacts,
        }
    finally:
        p.unlink(missing_ok=True)


# ── Plot helpers ───────────────────────────────────────────────────────────────

def _png_bfactor(bf_data, title):
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f: p = Path(f.name)
    plot_bfactor(bf_data, p, title=title)
    d = p.read_bytes(); p.unlink(missing_ok=True); return d

def _png_ramachandran(rama_data, title):
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f: p = Path(f.name)
    plot_ramachandran(rama_data, p, title=title)
    d = p.read_bytes(); p.unlink(missing_ok=True); return d

def _png_bsa_matrix(bsa_results, chain_ids, title):
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f: p = Path(f.name)
    plot_bsa_matrix(bsa_results, chain_ids, p, title=title)
    d = p.read_bytes(); p.unlink(missing_ok=True); return d

def _png_fingerprint(report):
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f: p = Path(f.name)
    plot_interaction_summary(report, p)
    d = p.read_bytes(); p.unlink(missing_ok=True); return d

def _png_2d(pdb_path, report, structure):
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f: p = Path(f.name)
    plot_ligand_2d(pdb_path, report, p, structure=structure)
    d = p.read_bytes(); p.unlink(missing_ok=True); return d


# ── CDR + Interface 3D Viewer ──────────────────────────────────────────────────

_EPITOPE_COLOR = "#00ccff"   # cyan  — antigen epitope residues

def _render_antibody_html(pdb_bytes: bytes, ab_data: dict,
                          width: int = 1100, height: int = 780) -> str:
    """
    3Dmol.js viewer with collapsible controls, sequence panel, and interface focus.
    Primary VH + VL + all Ag chains are shown; extra copies start collapsed.
    """
    pdb_str   = pdb_bytes.decode("utf-8", errors="replace")
    chains    = ab_data["chains"]
    vh_chains = ab_data["vh_chains"]
    vl_chains = ab_data["vl_chains"]
    ag_chains = ab_data["ag_chains"]
    pe        = ab_data.get("paratope_contacts", [])

    # Primary chains (first VH + first VL + all Ag); others start hidden
    vh_id = vh_chains[0] if vh_chains else None
    vl_id = vl_chains[0] if vl_chains else None
    primary = set(([vh_id] if vh_id else []) + ([vl_id] if vl_id else []) + ag_chains)
    all_chain_ids = list(chains.keys())

    # Chain type map & default visibility
    chain_types  = {cid: chains[cid]["type"] for cid in all_chain_ids}
    default_on   = {cid: (cid in primary) for cid in all_chain_ids}

    # ── CDR entries (all chains, but only primary shown by default) ────────────
    cdr_entries: list[dict] = []
    zoom_chain, zoom_resi = "", ""
    _zoom_map: dict = {}
    for cid, info in chains.items():
        for cdr_key, lbl_map in [
            ("cdr1", {"VH": "CDR-H1", "VL": "CDR-L1"}),
            ("cdr2", {"VH": "CDR-H2", "VL": "CDR-L2"}),
            ("cdr3", {"VH": "CDR-H3", "VL": "CDR-L3"}),
        ]:
            resis = info.get(f"{cdr_key}_resis", [])
            if resis and info["type"] in lbl_map:
                lbl = lbl_map[info["type"]]
                min_r, max_r = min(resis), max(resis)
                rng = f"{min_r}-{max_r}"
                cdr_entries.append({"chain": cid, "resi": rng,
                                    "label": lbl, "color": _CDR_COLORS[lbl]})
                if cid in primary:
                    _zoom_map[lbl] = (cid, rng)
    for lbl in ["CDR-H3", "CDR-L3", "CDR-H1", "CDR-L1", "CDR-H2", "CDR-L2"]:
        if lbl in _zoom_map:
            zoom_chain, zoom_resi = _zoom_map[lbl]; break

    # ── Epitope residues ───────────────────────────────────────────────────────
    epitope_resi_set: dict[str, set] = {}   # cid → set of pdb_resi
    for c in pe:
        epitope_resi_set.setdefault(c["Ag chain"], set()).add(c["Ag resi"])
    epitope_entries = [
        {"chain": cid, "resi": f"{min(rs)}-{max(rs)}"}
        for cid, rs in epitope_resi_set.items()
    ]

    # ── Sequence data for the sequence panel (primary chains only) ─────────────
    seq_rows: list[dict] = []
    for cid in ([vh_id] if vh_id else []) + ([vl_id] if vl_id else []) + ag_chains:
        if cid not in chains:
            continue
        info = chains[cid]
        ctype = info["type"]
        if ctype == "VH":
            cdr_resi_map = {"CDR-H1": set(info.get("cdr1_resis", [])),
                            "CDR-H2": set(info.get("cdr2_resis", [])),
                            "CDR-H3": set(info.get("cdr3_resis", []))}
        elif ctype == "VL":
            cdr_resi_map = {"CDR-L1": set(info.get("cdr1_resis", [])),
                            "CDR-L2": set(info.get("cdr2_resis", [])),
                            "CDR-L3": set(info.get("cdr3_resis", []))}
        else:
            cdr_resi_map = {}
        ep_resis = epitope_resi_set.get(cid, set())
        residues_js = []
        for pdb_resi, resname3, aa1 in info["residues"]:
            cdr_lbl = next((l for l, rs in cdr_resi_map.items() if pdb_resi in rs), None)
            residues_js.append({"r": pdb_resi, "a": aa1,
                                 "c": cdr_lbl, "e": pdb_resi in ep_resis})
        seq_rows.append({"id": cid, "type": ctype, "residues": residues_js})

    # ── JS data blobs ──────────────────────────────────────────────────────────
    cdr_entries_js     = json.dumps(cdr_entries)
    epitope_entries_js = json.dumps(epitope_entries)
    chain_types_js     = json.dumps(chain_types)
    default_on_js      = json.dumps(default_on)
    all_chain_ids_js   = json.dumps(all_chain_ids)
    seq_rows_js        = json.dumps(seq_rows)
    cdr_colors_js      = json.dumps(_CDR_COLORS)
    chain_colors_js    = json.dumps(_CHAIN_COLORS)
    zoom_chain_js      = json.dumps(zoom_chain)
    zoom_resi_js       = json.dumps(zoom_resi)
    epitope_color_js   = json.dumps(_EPITOPE_COLOR)

    pdb_escaped = pdb_str.replace("\\", "\\\\").replace("`", "\\`")

    viewer_h = height - 180   # leave room for controls + sequence

    return f"""<!DOCTYPE html>
<html><head>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.4/jquery.min.js"></script>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
* {{ box-sizing: border-box; }}
body {{ margin:0; background:#1a1a2e; color:#e8e8f0; font-family:sans-serif; font-size:12px; }}
#ab-wrap {{ display:flex; flex-direction:column; width:{width}px; }}

/* ── Controls ── */
.ab-ctrl-header {{
  display:flex; justify-content:space-between; align-items:center;
  background:#252540; padding:5px 10px; cursor:pointer;
  border-radius:4px; user-select:none;
}}
.ab-ctrl-header:hover {{ background:#2e2e55; }}
#ab-ctrl-body {{
  background:#1e1e38; padding:8px 12px;
  display:flex; flex-wrap:wrap; gap:12px; align-items:flex-start;
}}
.ab-ctrl-section {{ display:flex; flex-direction:column; gap:4px; min-width:120px; }}
.ab-ctrl-section b {{ font-size:11px; color:#aaa; text-transform:uppercase; letter-spacing:.05em; }}
.ab-ctrl-row {{ display:flex; flex-wrap:wrap; gap:6px; align-items:center; }}
label.ab-chk {{ display:flex; align-items:center; gap:4px; cursor:pointer; font-size:12px; }}
label.ab-chk input {{ cursor:pointer; }}
.ab-dim {{ opacity:0.45; }}

/* ── Sequence panel ── */
#ab-seq-panel {{ background:#1a1a2e; padding:4px 8px; max-height:110px; overflow-y:auto; }}
.ab-seq-row {{ display:flex; align-items:flex-start; gap:6px; margin-bottom:4px; }}
.ab-seq-lbl {{
  font-size:11px; font-weight:bold; white-space:nowrap;
  width:68px; flex-shrink:0; padding-top:2px;
}}
.ab-seq-body {{ display:flex; flex-wrap:wrap; gap:1px; }}
.ab-seq-aa {{
  display:inline-block; width:13px; text-align:center;
  font-size:10px; font-family:monospace; line-height:14px;
  border-radius:2px; cursor:pointer;
}}
.ab-seq-aa:hover {{ outline:1px solid #fff; }}
.ab-seq-sel {{ outline:2px solid #fff !important; }}

/* ── Viewer ── */
#ab-viewer {{ width:{width}px; height:{viewer_h}px; position:relative; }}

/* ── Surface opacity row ── */
#ab-surf-row {{ display:none; align-items:center; gap:6px; }}
</style>
</head><body>
<div id="ab-wrap">

<!-- Controls -->
<div class="ab-ctrl-header" onclick="abToggleControls()">
  <span>⚙ Controls</span><span id="ab-ctrl-arrow">▲</span>
</div>
<div id="ab-ctrl-body">

  <div class="ab-ctrl-section">
    <b>Chains</b>
    <div class="ab-ctrl-row" id="ab-chain-checks"></div>
  </div>

  <div class="ab-ctrl-section">
    <b>Style</b>
    <div class="ab-ctrl-row">
      <label class="ab-chk"><input type="radio" name="ab-style" value="cartoon" checked onchange="abRender()"> Cartoon</label>
      <label class="ab-chk"><input type="radio" name="ab-style" value="surface" onchange="abRender()"> Surface</label>
    </div>
    <div id="ab-surf-row">
      <span style="font-size:11px">Opacity</span>
      <input type="range" id="ab-surf-opacity" min="0.1" max="1.0" step="0.05" value="0.7"
             style="width:90px" oninput="abRender()">
      <span id="ab-surf-val">0.7</span>
    </div>
  </div>

  <div class="ab-ctrl-section">
    <b>CDRs</b>
    <div class="ab-ctrl-row" id="ab-cdr-checks"></div>
  </div>

</div><!-- end ctrl-body -->

<!-- Sequence panel -->
<div id="ab-seq-panel"></div>

<!-- 3D Viewer -->
<div id="ab-viewer"></div>
</div><!-- end wrap -->

<script>
$(function() {{
  // ── Data ─────────────────────────────────────────────────────────────────
  var ALL_CHAINS    = {all_chain_ids_js};
  var CHAIN_TYPES   = {chain_types_js};
  var DEFAULT_ON    = {default_on_js};
  var CDR_DATA      = {cdr_entries_js};
  var EPITOPE_DATA  = {epitope_entries_js};
  var SEQ_ROWS      = {seq_rows_js};
  var CDR_COLORS    = {cdr_colors_js};
  var CHAIN_COLORS  = {chain_colors_js};
  var EPITOPE_COLOR = {epitope_color_js};
  var ZOOM_CHAIN    = {zoom_chain_js};
  var ZOOM_RESI     = {zoom_resi_js};

  var surfObj = null;   // current surface object

  // ── Viewer init ───────────────────────────────────────────────────────────
  var viewer = $3Dmol.createViewer($("#ab-viewer"), {{backgroundColor:"0x1a1a2e"}});
  viewer.addModel(`{pdb_escaped}`, "pdb");

  // ── Build chain checkboxes ────────────────────────────────────────────────
  var chainChecks = document.getElementById("ab-chain-checks");
  ALL_CHAINS.forEach(function(cid) {{
    var ctype  = CHAIN_TYPES[cid] || "Unknown";
    var col    = CHAIN_COLORS[ctype] || "#888";
    var isPrim = DEFAULT_ON[cid];
    var lbl    = document.createElement("label");
    lbl.className = "ab-chk" + (isPrim ? "" : " ab-dim");
    lbl.innerHTML = '<input type="checkbox" id="ab-cc-' + cid + '"'
                  + (isPrim ? " checked" : "") + ' onchange="abRender()">'
                  + '<span style="color:' + col + '">' + cid
                  + ' <small>(' + ctype + ')</small></span>';
    chainChecks.appendChild(lbl);
  }});

  // ── Build CDR checkboxes ──────────────────────────────────────────────────
  var cdrChecks = document.getElementById("ab-cdr-checks");
  Object.keys(CDR_COLORS).forEach(function(lbl) {{
    var col = CDR_COLORS[lbl];
    var id  = "ab-cdr-" + lbl.replace(/-/g,"");
    var el  = document.createElement("label");
    el.className = "ab-chk";
    el.innerHTML = '<input type="checkbox" id="' + id + '" checked onchange="abRender()">'
                 + '<span style="color:' + col + '">' + lbl + '</span>';
    cdrChecks.appendChild(el);
  }});

  // ── Surface opacity live label ────────────────────────────────────────────
  document.getElementById("ab-surf-opacity").addEventListener("input", function() {{
    document.getElementById("ab-surf-val").textContent =
      parseFloat(this.value).toFixed(2);
  }});

  // ── Main render function ──────────────────────────────────────────────────
  window.abRender = function() {{
    var style = document.querySelector('input[name="ab-style"]:checked').value;
    var surfOpacity = parseFloat(document.getElementById("ab-surf-opacity").value);
    document.getElementById("ab-surf-row").style.display =
      style === "surface" ? "flex" : "none";

    // Remove previous surface
    if (surfObj !== null) {{ viewer.removeSurface(surfObj); surfObj = null; }}

    // Clear all styles
    viewer.setStyle({{}}, {{}});

    ALL_CHAINS.forEach(function(cid) {{
      var chk = document.getElementById("ab-cc-" + cid);
      if (!chk || !chk.checked) return;   // chain hidden

      var ctype   = CHAIN_TYPES[cid] || "Unknown";
      var col     = CHAIN_COLORS[ctype] || "#888888";
      var isPrim  = DEFAULT_ON[cid];
      var opacity = (ctype === "Antigen") ? 0.30 : 0.40;

      if (style === "cartoon") {{
        viewer.setStyle({{chain:cid}}, {{cartoon:{{color:col, opacity:opacity}}}});
      }} else {{
        // Surface: use base color, will add surface below
        viewer.setStyle({{chain:cid}}, {{cartoon:{{color:col, opacity:0.15}}}});
      }}
    }});

    // CDR loops
    CDR_DATA.forEach(function(e) {{
      var chainChk = document.getElementById("ab-cc-" + e.chain);
      var cdrChk   = document.getElementById("ab-cdr-" + e.label.replace(/-/g,""));
      if (chainChk && !chainChk.checked) return;
      if (cdrChk && !cdrChk.checked) return;
      viewer.setStyle({{chain:e.chain, resi:e.resi}},
                      {{cartoon:{{color:e.color, opacity:1.0}}}});
      viewer.addStyle({{chain:e.chain, resi:e.resi}},
                      {{stick:{{radius:0.22, colorscheme:"default"}}}});
    }});

    // Epitope residues
    EPITOPE_DATA.forEach(function(e) {{
      var chainChk = document.getElementById("ab-cc-" + e.chain);
      if (chainChk && !chainChk.checked) return;
      viewer.setStyle({{chain:e.chain, resi:e.resi}},
                      {{cartoon:{{color:EPITOPE_COLOR, opacity:1.0}}}});
      viewer.addStyle({{chain:e.chain, resi:e.resi}},
                      {{stick:{{radius:0.22, colorscheme:"default"}}}});
    }});

    // Surface
    if (style === "surface") {{
      var visChains = ALL_CHAINS.filter(function(c) {{
        var chk = document.getElementById("ab-cc-" + c); return chk && chk.checked;
      }});
      if (visChains.length > 0) {{
        surfObj = viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity:surfOpacity}},
          {{chain: visChains}});
      }}
    }}

    viewer.render();
    abColorSeq();
  }};

  // ── Sequence panel ────────────────────────────────────────────────────────
  function abBuildSeq() {{
    var panel = document.getElementById("ab-seq-panel");
    panel.innerHTML = "";
    SEQ_ROWS.forEach(function(row) {{
      var ctype = row.type;
      var col   = CHAIN_COLORS[ctype] || "#888";
      var div   = document.createElement("div");
      div.className = "ab-seq-row";
      var lbl = document.createElement("span");
      lbl.className = "ab-seq-lbl";
      lbl.style.color = col;
      lbl.textContent = ctype + " (" + row.id + ")";
      div.appendChild(lbl);
      var body = document.createElement("div");
      body.className = "ab-seq-body";
      row.residues.forEach(function(res) {{
        var span = document.createElement("span");
        span.className = "ab-seq-aa";
        span.textContent = res.a;
        span.dataset.chain = row.id;
        span.dataset.resi  = res.r;
        span.dataset.cdr   = res.c || "";
        span.dataset.ep    = res.e ? "1" : "";
        span.title = res.a + res.r + (res.c ? " [" + res.c + "]" : "") + (res.e ? " [epitope]" : "");
        span.onclick = function() {{
          document.querySelectorAll(".ab-seq-aa").forEach(function(s) {{ s.classList.remove("ab-seq-sel"); }});
          span.classList.add("ab-seq-sel");
          viewer.zoomTo({{chain:span.dataset.chain, resi:parseInt(span.dataset.resi)}});
          viewer.render();
        }};
        body.appendChild(span);
      }});
      div.appendChild(body);
      panel.appendChild(div);
    }});
    abColorSeq();
  }}

  window.abColorSeq = function() {{
    document.querySelectorAll(".ab-seq-aa").forEach(function(span) {{
      if (span.classList.contains("ab-seq-sel")) return;
      var cdr = span.dataset.cdr;
      var ep  = span.dataset.ep;
      if (cdr && CDR_COLORS[cdr]) {{
        var c = CDR_COLORS[cdr];
        span.style.background = c + "44";
        span.style.color = c;
      }} else if (ep) {{
        span.style.background = EPITOPE_COLOR + "44";
        span.style.color = EPITOPE_COLOR;
      }} else {{
        span.style.background = "";
        span.style.color = "#ccc";
      }}
    }});
  }};

  // ── Collapsible controls ──────────────────────────────────────────────────
  window.abToggleControls = function() {{
    var body  = document.getElementById("ab-ctrl-body");
    var arrow = document.getElementById("ab-ctrl-arrow");
    if (body.style.display === "none") {{
      body.style.display = "flex"; arrow.textContent = "▲";
    }} else {{
      body.style.display = "none"; arrow.textContent = "▼";
    }}
  }};

  // ── Initial render ────────────────────────────────────────────────────────
  abBuildSeq();
  abRender();
  if (ZOOM_CHAIN && ZOOM_RESI) {{
    viewer.zoomTo({{chain:ZOOM_CHAIN, resi:ZOOM_RESI}});
  }} else {{
    viewer.zoomTo();
  }}
  viewer.render();
}});
</script>
</body></html>"""


# ═══════════════════════════════════════════════════════════════════════════════
# File upload + top-level controls
# ═══════════════════════════════════════════════════════════════════════════════

up_col, opt_col = st.columns([11, 1])
with up_col:
    uploaded = st.file_uploader(
        "Upload a PDB file",
        type=["pdb", "ent"],
        help="Standard PDB format (.pdb or .ent)",
        label_visibility="collapsed",
    )
with opt_col:
    with st.popover("⚙", help="Analysis settings"):
        contact_cutoff = st.slider(
            "Contact cutoff (Å)", min_value=3.0, max_value=7.0,
            value=4.5, step=0.1,
        )

st.divider()

if uploaded:
    pdb_bytes = uploaded.read()

    if st.session_state.get("pdb_name") != uploaded.name:
        d = Path(tempfile.mkdtemp())
        p = d / uploaded.name
        p.write_bytes(pdb_bytes)
        st.session_state["pdb_name"] = uploaded.name
        st.session_state["pdb_path"] = str(p)

    pdb_path = Path(st.session_state["pdb_path"])
    st.title(f"🧫 {pdb_path.stem.upper()}")
    st.caption(f"File: `{uploaded.name}`  ·  {len(pdb_bytes):,} bytes")

    with st.spinner("Running structure analysis…"):
        struct_data = _run_structure_analysis(pdb_bytes)

    summary      = struct_data["summary"]
    bf_data      = struct_data["bfactor"]
    rama_data    = struct_data["ramachandran"]
    rg           = struct_data["rg"]
    fasta_str    = struct_data["fasta"]
    missing_data = struct_data["missing"]

    with st.spinner("Running ligand interaction analysis…"):
        reports = _run_ligand_analysis(pdb_bytes, contact_cutoff)

    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Resolution", f"{summary['resolution']} Å" if summary['resolution'] else "N/A")
    c2.metric("R-work",     str(summary["r_work"])  if summary["r_work"]  else "N/A")
    c3.metric("R-free",     str(summary["r_free"])  if summary["r_free"]  else "N/A")
    c4.metric("Residues",   f"{summary['n_residues']:,}")
    c5.metric("Ligands",    summary["n_ligands"])
    c6.metric("Rg (Å)",     f"{rg:.2f}" if rg else "N/A")
    st.divider()

else:
    pdb_bytes = pdb_path = summary = bf_data = rama_data = None
    rg = fasta_str = missing_data = reports = None


# ═══════════════════════════════════════════════════════════════════════════════
# Tabs
# ═══════════════════════════════════════════════════════════════════════════════
tab_summary, tab_bsa, tab_viewer, tab_2d, tab_ab = st.tabs([
    "📋 Structure Summary",
    "⬛ Buried Surface Area",
    "🌐 3D Viewer",
    "🖼 2D Diagrams",
    "🧫 Antibody Analysis",
])


# ─── TAB 1: Structure Summary ──────────────────────────────────────────────────
with tab_summary:
    if not uploaded:
        st.info(_NO_PDB)
    else:
        st.subheader("Metadata")
        meta = {
            "PDB ID":             summary["pdb_id"],
            "Title":              summary["title"] or "—",
            "Method":             summary["method"] or "—",
            "Organism":           summary["organism"] or "—",
            "Deposited":          summary["deposition_date"] or "—",
            "Models":             summary["n_models"],
            "Chains":             summary["chains"],
            "Residues per chain": summary["residues_per_chain"],
            "Atoms":              f"{summary['n_atoms']:,}",
            "Ligand names":       summary["ligand_names"] or "—",
            "Waters":             summary["n_waters"],
        }
        st.dataframe(
            pd.DataFrame(meta.items(), columns=["Field", "Value"]),
            use_container_width=True, hide_index=True,
        )

        if missing_data:
            st.warning(f"**Missing residue gaps** detected in {len(missing_data)} chain(s)")
            st.dataframe(pd.DataFrame(missing_data), use_container_width=True, hide_index=True)
        else:
            st.success("No missing residue gaps detected")

        all_ligs = _list_all_ligands(pdb_bytes)
        if all_ligs:
            with st.expander(f"All HETATM residues ({len(all_ligs)})"):
                st.dataframe(pd.DataFrame(all_ligs), use_container_width=True, hide_index=True)

        st.divider()
        st.subheader("B-factor Profile")
        if bf_data:
            n_high = sum(1 for r in bf_data if r["flag"] == "high")
            n_low  = sum(1 for r in bf_data if r["flag"] == "low")
            st.caption(f"{len(bf_data)} Cα residues  ·  **{n_high}** high (>2σ)  ·  **{n_low}** low (<2σ)")
            bf_png = _png_bfactor(bf_data, pdb_path.stem.upper())
            st.image(bf_png, width="stretch")
            col_dl, _ = st.columns([1, 5])
            with col_dl:
                st.download_button("⬇ B-factor PNG", bf_png,
                                   file_name=f"{pdb_path.stem}_bfactor.png", mime="image/png")
            with st.expander("B-factor data table"):
                st.dataframe(pd.DataFrame(bf_data), use_container_width=True, hide_index=True)
        else:
            st.info("No B-factor data (no Cα atoms found).")

        st.divider()
        st.subheader("Ramachandran Plot")
        if rama_data:
            n_total   = len(rama_data)
            n_outlier = sum(1 for r in rama_data if r["region"] == "outlier")
            n_helix   = sum(1 for r in rama_data if r["region"] == "helix")
            n_sheet   = sum(1 for r in rama_data if r["region"] == "sheet")
            pct_ok    = 100 * (1 - n_outlier / n_total) if n_total else 0
            rc1, rc2, rc3, rc4 = st.columns(4)
            rc1.metric("Total residues", n_total)
            rc2.metric("Helix", n_helix)
            rc3.metric("Sheet", n_sheet)
            rc4.metric("Outliers", n_outlier, delta=f"{pct_ok:.1f}% OK", delta_color="normal")
            rama_png = _png_ramachandran(rama_data, pdb_path.stem.upper())
            col_img, col_dl = st.columns([4, 1])
            with col_img:
                st.image(rama_png)
            with col_dl:
                st.download_button("⬇ Ramachandran PNG", rama_png,
                                   file_name=f"{pdb_path.stem}_rama.png", mime="image/png")
            if n_outlier:
                with st.expander(f"Outlier residues ({n_outlier})"):
                    st.dataframe(
                        pd.DataFrame([r for r in rama_data if r["region"] == "outlier"]),
                        use_container_width=True, hide_index=True,
                    )
        else:
            st.info("No Ramachandran data computed.")

        st.divider()
        st.subheader("FASTA Sequences (from ATOM records)")
        if fasta_str:
            st.code(fasta_str, language=None)
            st.download_button("⬇ Download FASTA", fasta_str,
                               file_name=f"{pdb_path.stem}.fasta", mime="text/plain")
        else:
            st.info("No protein sequence extracted.")


# ─── TAB 2: Buried Surface Area ───────────────────────────────────────────────
with tab_bsa:
    if not uploaded:
        st.info(_NO_PDB)
    else:
        st.subheader("Buried Surface Area — Chain Interfaces")
        _, bsa_opt_col = st.columns([11, 1])
        with bsa_opt_col:
            with st.popover("⚙", help="BSA settings"):
                probe_radius = st.slider(
                    "Probe radius (Å)", min_value=1.0, max_value=2.0,
                    value=1.4, step=0.1,
                )
                n_points_bsa = st.select_slider(
                    "Sphere points", options=[50, 100, 200, 500], value=100,
                )
        st.caption("BSA(A, B) = ( SASA_A + SASA_B − SASA_AB ) / 2")
        n_chains = len(summary["chains"].split(",")) if summary["chains"] else 0
        if n_chains < 2:
            st.warning("Only one chain detected — BSA requires at least two chains.")
        else:
            with st.spinner(f"Computing BSA (probe={probe_radius} Å)…"):
                bsa_results, chain_ids = _run_bsa(pdb_bytes, probe_radius, n_points_bsa)
            if not bsa_results:
                st.info("No inter-chain contacts found.")
            else:
                top = bsa_results[0]
                n_iface = sum(1 for r in bsa_results if r["bsa_total"] > 0)
                st.success(
                    f"**{n_iface}** interface(s)  ·  Largest: chains "
                    f"{top['chain_1']}–{top['chain_2']} ({top['bsa_total']:.0f} Å²)"
                )
                bsa_png = _png_bsa_matrix(bsa_results, chain_ids, pdb_path.stem.upper())
                st.image(bsa_png, width="stretch")
                st.download_button("⬇ BSA matrix PNG", bsa_png,
                                   file_name=f"{pdb_path.stem}_bsa_matrix.png", mime="image/png")
                st.divider()
                df_bsa = pd.DataFrame(bsa_results).rename(columns={
                    "chain_1": "Chain 1", "chain_2": "Chain 2",
                    "sasa_1": "SASA 1 (A\u00b2)", "sasa_2": "SASA 2 (A\u00b2)",
                    "sasa_complex": "SASA complex (A\u00b2)",
                    "bsa_total": "BSA total (A\u00b2)",
                    "bsa_on_1": "BSA on 1 (A\u00b2)", "bsa_on_2": "BSA on 2 (A\u00b2)",
                    "interface_pct_1": "Interface % 1", "interface_pct_2": "Interface % 2",
                })
                bsa_col = "BSA total (A\u00b2)"
                st.dataframe(
                    df_bsa.style.background_gradient(subset=[bsa_col], cmap="YlOrRd"),
                    hide_index=True,
                )
                st.download_button("⬇ BSA table CSV", df_bsa.to_csv(index=False),
                                   file_name=f"{pdb_path.stem}_bsa.csv", mime="text/csv")


# ─── TAB 3: 3D Viewer + Ligand Analysis ───────────────────────────────────────
with tab_viewer:
    if not uploaded:
        st.info(_NO_PDB)
    else:
        st.subheader("Interactive 3D Viewer")
        if not reports:
            # No small-molecule ligand — show the antibody interface viewer instead
            with st.spinner("Building antibody interface viewer…"):
                ab_data_3d = _run_antibody_analysis(pdb_bytes)
            _ab_chains    = ab_data_3d["vh_chains"] + ab_data_3d["vl_chains"]
            _ag_chains_3d = ab_data_3d["ag_chains"]
            if _ab_chains:
                _id_parts = []
                if ab_data_3d["vh_chains"]:
                    _id_parts.append(f"VH: {', '.join(ab_data_3d['vh_chains'])}")
                if ab_data_3d["vl_chains"]:
                    _id_parts.append(f"VL: {', '.join(ab_data_3d['vl_chains'])}")
                if _ag_chains_3d:
                    _id_parts.append(f"Antigen: {', '.join(_ag_chains_3d)}")
                st.caption("  ·  ".join(_id_parts) + "  ·  CDR loops coloured · gold = epitope · drag to rotate · scroll to zoom")
                _ab_html = _render_antibody_html(pdb_bytes, ab_data_3d, width=1100, height=780)
                components.html(_ab_html, height=820, scrolling=False)
                st.download_button(
                    "⬇ Download interface viewer HTML", _ab_html,
                    file_name=f"{pdb_path.stem}_antibody.html", mime="text/html",
                )
            else:
                st.info("No antibody chains or ligands detected in this structure.")
        else:
            lig_labels = [f"{r.ligand_resn}  chain {r.ligand_chain}:{r.ligand_resi}" for r in reports]
            sel_idx = st.selectbox("Select ligand", range(len(reports)),
                                   format_func=lambda i: lig_labels[i]) if len(reports) > 1 else 0
            if len(reports) == 1:
                st.caption(f"Showing: **{lig_labels[0]}**")
            report = reports[sel_idx]
            st.caption(f"{report.summary()}  ·  Drag to rotate · scroll to zoom.")
            lc1, lc2, lc3, lc4, lc5 = st.columns(5)
            lc1.markdown("🟡 **Ligand**"); lc2.markdown("🔵 **Contacts**")
            lc3.markdown("🟠 **H-bonds**"); lc4.markdown("🟣 **Pi**")
            lc5.markdown("🔵🔴 **Salt bridges**")

            @st.cache_resource
            def _cached_structure(path_str: str):
                return al_load(path_str)

            structure = _cached_structure(str(pdb_path))
            html_str  = _render_html(pdb_path, report, structure, width=1100, height=680)
            components.html(html_str, height=720, scrolling=False)
            st.download_button(
                "⬇ Download standalone HTML viewer", html_str,
                file_name=f"{report.ligand_resn}_{report.ligand_chain}{report.ligand_resi}.html",
                mime="text/html",
            )

        st.divider()
        st.subheader("Ligand Interaction Analysis")
        if not reports:
            st.warning("No primary ligands detected.")
        else:
            for rep in reports:
                lig_id = f"{rep.ligand_resn} chain {rep.ligand_chain} resi {rep.ligand_resi}"
                counts = (f"{rep.n_contacts} contacts · {rep.n_hbonds} H-bonds · "
                          f"{rep.n_pi} pi · {rep.n_salt_bridges} salt bridges")
                with st.expander(f"**{lig_id}**  —  {counts}", expanded=False):
                    it1, it2, it3, it4 = st.tabs([
                        f"Contacts ({rep.n_contacts})", f"H-bonds ({rep.n_hbonds})",
                        f"Pi ({rep.n_pi})", f"Salt ({rep.n_salt_bridges})",
                    ])
                    key = f"{rep.ligand_resn}_{rep.ligand_resi}"
                    with it1:
                        if rep.contacts:
                            df = pd.DataFrame([c.as_dict() for c in rep.contacts])
                            st.dataframe(df, use_container_width=True, hide_index=True)
                            st.download_button("⬇ CSV", df.to_csv(index=False),
                                               file_name=f"{rep.ligand_resn}_contacts.csv",
                                               mime="text/csv", key=f"c_{key}")
                        else:
                            st.info("No contacts detected.")
                    with it2:
                        if rep.hbonds:
                            df = pd.DataFrame([h.as_dict() for h in rep.hbonds])
                            st.dataframe(df, use_container_width=True, hide_index=True)
                        else:
                            st.info("No hydrogen bonds detected.")
                    with it3:
                        if rep.pi_interactions:
                            df = pd.DataFrame([p_.as_dict() for p_ in rep.pi_interactions])
                            st.dataframe(df, use_container_width=True, hide_index=True)
                        else:
                            st.info("No pi interactions detected.")
                    with it4:
                        if rep.salt_bridges:
                            df = pd.DataFrame([s.as_dict() for s in rep.salt_bridges])
                            st.dataframe(df, use_container_width=True, hide_index=True)
                        else:
                            st.info("No salt bridges detected.")


# ─── TAB 4: 2D Diagrams ───────────────────────────────────────────────────────
with tab_2d:
    if not uploaded:
        st.info(_NO_PDB)
    else:
        st.subheader("2D Interaction Diagrams")
        if not reports:
            st.warning("No primary ligands detected.")
        else:
            @st.cache_resource
            def _cached_structure_2d(path_str: str):
                return al_load(path_str)

            structure_2d = _cached_structure_2d(str(pdb_path))
            for rep in reports:
                lig_label = f"{rep.ligand_resn} chain {rep.ligand_chain}:{rep.ligand_resi}"
                st.markdown(f"### {lig_label}")
                col_fp, col_2d = st.columns(2, gap="large")
                key = f"{rep.ligand_resn}_{rep.ligand_resi}"
                has_int = (rep.n_contacts + rep.n_hbonds + rep.n_pi + rep.n_salt_bridges) > 0
                with col_fp:
                    st.markdown("#### Interaction Fingerprint")
                    if has_int:
                        fp = _png_fingerprint(rep)
                        st.image(fp, width="stretch")
                        st.download_button("⬇ Fingerprint PNG", fp,
                                           file_name=f"{rep.ligand_resn}_fingerprint.png",
                                           mime="image/png", key=f"fp_{key}")
                    else:
                        st.info("No interactions to plot.")
                with col_2d:
                    st.markdown("#### LIGPLOT-style 2D Diagram")
                    if has_int:
                        d2 = _png_2d(pdb_path, rep, structure_2d)
                        st.image(d2, width="stretch")
                        st.download_button("⬇ 2D diagram PNG", d2,
                                           file_name=f"{rep.ligand_resn}_2d.png",
                                           mime="image/png", key=f"2d_{key}")
                    else:
                        st.info("No interactions to plot.")
                st.divider()


# ─── TAB 5: Antibody Analysis ─────────────────────────────────────────────────
with tab_ab:
    if not uploaded:
        st.info(_NO_PDB)
    else:
        with st.spinner("Running antibody chain analysis…"):
            ab_data = _run_antibody_analysis(pdb_bytes)

        chains    = ab_data["chains"]
        vh_chains = ab_data["vh_chains"]
        vl_chains = ab_data["vl_chains"]
        ag_chains = ab_data["ag_chains"]

        if not chains:
            st.warning("No protein chains with ≥30 residues found.")
        elif not vh_chains and not vl_chains:
            st.warning(
                "No antibody variable domains detected.  \n"
                "The heuristic looks for conserved FR2 tryptophan motifs "
                "(W-[VIL]-[RK] for VH, W-[YF]-QQ for VL). "
                "This may not work for very unusual sequences."
            )
        else:
            # ── Chain classification summary ───────────────────────────────────
            st.subheader("Chain Classification")
            st.caption(
                "Heuristic classification — no HMMER required. "
                "CDR boundaries are approximate (±2 residues, Chothia/Kabat offsets from conserved anchors)."
            )

            am1, am2, am3, am4 = st.columns(4)
            am1.metric("VH chains",      len(vh_chains))
            am2.metric("VL chains",      len(vl_chains))
            am3.metric("Antigen chains", len(ag_chains))
            cdrs_found = sum(
                1 for cid in chains
                for k in ("cdr1_seq", "cdr2_seq", "cdr3_seq")
                if chains[cid].get(k)
            )
            am4.metric("CDRs detected", cdrs_found)

            chain_rows = []
            for cid, info in chains.items():
                chain_rows.append({
                    "Chain": cid,
                    "Type":  info["type"],
                    "Length (AA)": info["length"],
                    "CDR1": info.get("cdr1_seq") or "—",
                    "CDR2": info.get("cdr2_seq") or "—",
                    "CDR3": info.get("cdr3_seq") or "—",
                })
            df_chains = pd.DataFrame(chain_rows)

            def _color_type(val: str) -> str:
                colors = {"VH": "background-color:#1a3a5c;color:#7eb8e8",
                          "VL": "background-color:#5c3a00;color:#f0b060",
                          "Unknown": ""}
                return colors.get(val, "")

            st.dataframe(
                df_chains.style.applymap(_color_type, subset=["Type"]),
                use_container_width=True, hide_index=True,
            )

            # ── CDR detail table ───────────────────────────────────────────────
            st.divider()
            st.subheader("CDR Sequences & Lengths")
            cdr_rows = []
            for cid in vh_chains + vl_chains:
                info   = chains[cid]
                prefix = "H" if info["type"] == "VH" else "L"
                for n_, lbl in [(1, f"CDR-{prefix}1"), (2, f"CDR-{prefix}2"), (3, f"CDR-{prefix}3")]:
                    seq_  = info.get(f"cdr{n_}_seq")
                    resis = info.get(f"cdr{n_}_resis", [])
                    cdr_rows.append({
                        "Chain": cid,
                        "Type":  info["type"],
                        "CDR":   lbl,
                        "Sequence": seq_ or "—",
                        "Length": len(seq_) if seq_ else 0,
                        "Resi range": f"{min(resis)}–{max(resis)}" if resis else "—",
                    })
            if cdr_rows:
                df_cdr = pd.DataFrame(cdr_rows)
                st.dataframe(df_cdr, use_container_width=True, hide_index=True)
                st.download_button(
                    "⬇ CDR table CSV", df_cdr.to_csv(index=False),
                    file_name=f"{pdb_path.stem}_cdrs.csv", mime="text/csv",
                )

            # ── Interface 3D Viewer ────────────────────────────────────────────
            st.divider()
            st.subheader("Interface 3D Viewer")

            # Summarise what was auto-detected
            _id_parts = []
            if vh_chains:
                _id_parts.append(f"**VH** chain(s): {', '.join(vh_chains)}")
            if vl_chains:
                _id_parts.append(f"**VL** chain(s): {', '.join(vl_chains)}")
            if ag_chains:
                _id_parts.append(f"**Antigen** chain(s): {', '.join(ag_chains)}")
            else:
                _id_parts.append("no antigen chain detected")
            st.caption("  ·  ".join(_id_parts))

            _has_pe = bool(ab_data.get("paratope_contacts"))
            if _has_pe:
                _n_ep = len({(c["Ag chain"], c["Ag resi"])
                              for c in ab_data["paratope_contacts"]})
                st.caption(
                    f"View centred on the paratope tip (CDR-H3 preferred).  "
                    f"Gold = **{_n_ep}** epitope residues.  "
                    f"Framework is semi-transparent to expose the binding site."
                )
            else:
                st.caption(
                    "View centred on CDR loops.  "
                    "No antigen contacts to highlight — upload a co-crystal structure "
                    "to see the epitope."
                )

            ab_html = _render_antibody_html(pdb_bytes, ab_data, width=1100, height=780)
            components.html(ab_html, height=820, scrolling=False)
            st.download_button(
                "⬇ Download interface viewer HTML", ab_html,
                file_name=f"{pdb_path.stem}_antibody.html", mime="text/html",
            )

            # ── VH-VL Interface ────────────────────────────────────────────────
            if vh_chains and vl_chains:
                st.divider()
                st.subheader(f"VH–VL Interface  (chains {vh_chains[0]}–{vl_chains[0]})")
                vhvl = ab_data["vhvl_contacts"]
                if vhvl:
                    df_vhvl = pd.DataFrame(vhvl)
                    st.caption(f"{len(vhvl)} unique residue-pair contacts ≤ 4.5 Å")
                    dist_col = "dist (Å)"
                    st.dataframe(
                        df_vhvl.style.background_gradient(subset=[dist_col], cmap="Blues_r"),
                        hide_index=True,
                    )
                    st.download_button(
                        "⬇ VH-VL contacts CSV", df_vhvl.to_csv(index=False),
                        file_name=f"{pdb_path.stem}_vhvl_contacts.csv", mime="text/csv",
                    )
                else:
                    st.info("No VH-VL contacts found within 4.5 Å.")

            # ── Paratope-Epitope ───────────────────────────────────────────────
            if ag_chains:
                st.divider()
                st.subheader(f"Paratope–Epitope Contacts  (antigen: {', '.join(ag_chains)})")
                pe = ab_data["paratope_contacts"]
                if pe:
                    df_pe = pd.DataFrame(pe)
                    cdr_count = sum(1 for r in pe if r["In CDR"] == "✓")
                    st.caption(
                        f"{len(pe)} unique contacts ≤ 4.5 Å  ·  "
                        f"**{cdr_count}** in CDR residues"
                    )
                    dist_col = "dist (Å)"
                    st.dataframe(
                        df_pe.style.background_gradient(subset=[dist_col], cmap="Oranges_r"),
                        hide_index=True,
                    )
                    st.download_button(
                        "⬇ Paratope-epitope CSV", df_pe.to_csv(index=False),
                        file_name=f"{pdb_path.stem}_paratope_epitope.csv", mime="text/csv",
                    )
                else:
                    st.info("No antibody–antigen contacts found within 4.5 Å.")
            elif vh_chains or vl_chains:
                st.divider()
                st.info(
                    "No antigen chain detected. All protein chains were classified as VH or VL.  \n"
                    "If an antigen is present but undetected, it may have been classified as VH/VL "
                    "due to a coincidental FR2-like motif — check the chain classification table above."
                )
