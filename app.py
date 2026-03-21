#!/usr/bin/env python3
"""
app.py — Streamlit web app for PDB structure analysis.

Integrates:
  • summarize_structures.py  — metadata, B-factors, Ramachandran, FASTA, Rg
  • analyze_ligands.py       — contacts, H-bonds, pi interactions, salt bridges
  • visualize_interactions.py — 3D HTML viewer (3Dmol.js) + 2D LIGPLOT diagrams

Run:
    streamlit run app.py
"""

from __future__ import annotations

import json
import os
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

# ── Add toolkit directory to path ─────────────────────────────────────────────
TOOLKIT_DIR = Path(__file__).parent
_FILES_DIR = TOOLKIT_DIR / "Files"
_HUMSAVAR_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/"
    "current_release/knowledgebase/variants/humsavar.txt"
)
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

# ── Page config ────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="PDB Structure Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Custom CSS ─────────────────────────────────────────────────────────────────
st.markdown("""
<style>
.stTabs [data-baseweb="tab-list"] { gap: 6px; }
.stTabs [data-baseweb="tab"] { padding: 8px 18px; border-radius: 6px 6px 0 0; }
</style>
""", unsafe_allow_html=True)

_NO_PDB = "Upload a PDB file in the sidebar to get started."


# ═══════════════════════════════════════════════════════════════════════════════
# Cached analysis functions
# ═══════════════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def _run_structure_analysis(data: bytes) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data)
        p = Path(f.name)
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
        f.write(data)
        p = Path(f.name)
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
        return analyze_ligand(
            structure, " or ".join(parts),
            pdb_path=p, contact_cutoff=cutoff,
        )
    finally:
        p.unlink(missing_ok=True)


@st.cache_data(show_spinner=False)
def _run_bsa(data: bytes, probe_radius: float, n_points: int):
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data)
        p = Path(f.name)
    try:
        structure = ss_load(str(p))
        chain_ids = [ch.id for ch in structure[0]]
        results   = buried_surface_areas(structure,
                                         probe_radius=probe_radius,
                                         n_points=n_points)
        return results, chain_ids
    finally:
        p.unlink(missing_ok=True)


@st.cache_data(show_spinner=False)
def _list_all_ligands(data: bytes) -> list[dict]:
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data)
        p = Path(f.name)
    try:
        structure = al_load(str(p))
        ligs = get_ligands(structure, exclude_water=True)
        return [
            {
                "resn":  res.get_resname().strip(),
                "chain": res.get_parent().id,
                "resi":  res.id[1],
                "atoms": sum(1 for _ in res.get_atoms()),
            }
            for res in ligs
        ]
    finally:
        p.unlink(missing_ok=True)


# ═══════════════════════════════════════════════════════════════════════════════
# AlphaFold + UniProt helpers
# ═══════════════════════════════════════════════════════════════════════════════


@st.cache_data(show_spinner=False)
def _fetch_alphafold_pdb(uniprot_id: str) -> bytes | None:
    uid = uniprot_id.upper()
    local_path = _FILES_DIR / f"AF-{uid}.pdb"
    if local_path.exists():
        return local_path.read_bytes()
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
    try:
        req = urllib.request.Request(api_url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=20) as r:
            entries = json.loads(r.read())
        entry = entries[0] if isinstance(entries, list) else entries
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            return None
    except Exception:
        return None
    try:
        with urllib.request.urlopen(pdb_url, timeout=30) as r:
            data = r.read()
        _FILES_DIR.mkdir(exist_ok=True)
        local_path.write_bytes(data)
        return data
    except Exception:
        return None


@st.cache_data(show_spinner=False)
def _fetch_uniprot(uniprot_id: str) -> dict | None:
    """Fetch UniProt JSON for structural features (disulfide bonds, glycosylation, etc.).
    Downloads once to Files/{uid}_uniprot.json and reads locally on subsequent calls."""
    uid = uniprot_id.upper()
    local_path = _FILES_DIR / f"{uid}_uniprot.json"
    if local_path.exists():
        return json.loads(local_path.read_text())
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=20) as r:
            raw = r.read()
        _FILES_DIR.mkdir(exist_ok=True)
        local_path.write_bytes(raw)
        return json.loads(raw)
    except Exception:
        return None


@st.cache_data(show_spinner=False)
def _load_humsavar() -> dict[str, list[dict]]:
    """Download humsavar.txt once to Files/ and parse into a dict keyed by UniProt AC.

    Each value is a list of variant dicts with keys:
      position, mutation, description, disease, is_disease, dbsnp, category
    """
    local_path = _FILES_DIR / "humsavar.txt"
    if not local_path.exists():
        _FILES_DIR.mkdir(exist_ok=True)
        with urllib.request.urlopen(_HUMSAVAR_URL, timeout=120) as r:
            local_path.write_bytes(r.read())

    variants: dict[str, list[dict]] = {}
    with local_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            # Skip blank lines, comment/header lines, and the column-header row
            if not line or line.startswith("//") or line.startswith("Gene"):
                continue
            parts = line.split(None, 6)
            if len(parts) < 6:
                continue
            gene, ac, ftid, aa_change, category, dbsnp = parts[:6]
            disease_name = parts[6].strip() if len(parts) > 6 else ""
            if disease_name == "-":
                disease_name = ""
            m = re.search(r"\d+", aa_change)
            position = int(m.group()) if m else None
            if position is None:
                continue
            is_disease = category == "LP/P"
            variants.setdefault(ac, []).append({
                "position":    position,
                "mutation":    aa_change,
                "description": disease_name or aa_change,
                "disease":     disease_name if is_disease else "",
                "is_disease":  is_disease,
                "dbsnp":       dbsnp if dbsnp != "-" else "",
                "clinvar":     "",
                "category":    category,
            })
    return variants



def _parse_uniprot_features(data: dict) -> dict:
    result: dict = {
        "disulfide_bonds": [], "glycosylation": [],
        "active_sites": [],    "binding_sites": [],
        "modified_residues": [], "signal_peptide": [], "transmembrane": [],
    }
    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        loc   = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end   = loc.get("end",   {}).get("value")
        desc  = feat.get("description", "")
        if start is None:
            continue
        if ftype == "Disulfide bond" and end is not None:
            result["disulfide_bonds"].append({"pos1": start, "pos2": end})
        elif ftype == "Glycosylation":
            result["glycosylation"].append({"position": start, "glycan_type": desc or "glycosylation"})
        elif ftype == "Active site":
            result["active_sites"].append({"position": start, "description": desc or "active site"})
        elif ftype == "Binding site":
            result["binding_sites"].append({"position": start, "description": desc or "binding site"})
        elif ftype == "Modified residue":
            result["modified_residues"].append({"position": start, "description": desc})
        elif ftype == "Signal peptide":
            result["signal_peptide"].append({"start": start, "end": end})
        elif ftype == "Transmembrane":
            result["transmembrane"].append({"start": start, "end": end})
    return result


# ─── Plot helpers ─────────────────────────────────────────────────────────────

def _png_bfactor(bf_data: list[dict], title: str) -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_bfactor(bf_data, p, title=title)
    data = p.read_bytes(); p.unlink(missing_ok=True); return data


def _png_ramachandran(rama_data: list[dict], title: str) -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_ramachandran(rama_data, p, title=title)
    data = p.read_bytes(); p.unlink(missing_ok=True); return data


def _png_bsa_matrix(bsa_results: list, chain_ids: list, title: str) -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_bsa_matrix(bsa_results, chain_ids, p, title=title)
    data = p.read_bytes(); p.unlink(missing_ok=True); return data


def _png_fingerprint(report, title: str = "") -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_interaction_summary(report, p)
    data = p.read_bytes(); p.unlink(missing_ok=True); return data


def _png_2d(pdb_path: Path, report, structure) -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_ligand_2d(pdb_path, report, p, structure=structure)
    data = p.read_bytes(); p.unlink(missing_ok=True); return data


# ═══════════════════════════════════════════════════════════════════════════════
# AlphaFold HTML renderer
# ═══════════════════════════════════════════════════════════════════════════════

_AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def _render_alphafold_html(pdb_bytes: bytes,
                           variants: list[dict],
                           features: dict | None = None,
                           protein_name: str = "",
                           width: int = 1100,
                           height: int = 680) -> str:
    if features is None:
        features = {}

    pdb_str = pdb_bytes.decode("utf-8", errors="replace").replace("`", "'")

    disease_vars    = [v for v in variants if v["is_disease"]]
    other_vars      = [v for v in variants if not v["is_disease"]]
    disulfide_bonds = features.get("disulfide_bonds", [])
    glycosylation   = features.get("glycosylation",   [])
    active_sites    = features.get("active_sites",    [])
    binding_sites   = features.get("binding_sites",   [])

    # Extract sequence from CA atoms (keyed by chain+resi to handle multi-chain)
    seq_residues: list[dict] = []
    seen_keys: set[tuple] = set()
    for line in pdb_bytes.decode("utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") and len(line) > 26 and line[12:16].strip() == "CA":
            try:
                chain = line[21]
                resi  = int(line[22:26].strip())
                resn  = line[17:20].strip()
                key   = (chain, resi)
                if key not in seen_keys:
                    seen_keys.add(key)
                    seq_residues.append({"resi": resi, "resn": resn,
                                         "aa": _AA3TO1.get(resn, "?"), "chain": chain})
            except (ValueError, IndexError):
                pass
    seq_js = json.dumps(seq_residues)

    def _var_js(vlist: list[dict]) -> str:
        return json.dumps([
            {"resi": v["position"], "label": v["mutation"],
             "disease": v["disease"], "desc": v["description"][:120]}
            for v in vlist if v.get("position")
        ])

    disease_js   = _var_js(disease_vars)
    other_js     = _var_js(other_vars)
    disulfide_js = json.dumps([{"pos1": b["pos1"], "pos2": b["pos2"]} for b in disulfide_bonds])
    glycan_js    = json.dumps([{"resi": g["position"], "desc": g["glycan_type"][:80]} for g in glycosylation])
    active_js    = json.dumps([{"resi": a["position"], "desc": a["description"][:80]} for a in active_sites])
    binding_js   = json.dumps([{"resi": b["position"], "desc": b["description"][:80]} for b in binding_sites])

    def _chk(id_: str, color: str, label: str, count: int) -> str:
        if count == 0:
            return ""
        return (f'<label class="ctrl-row"><input type="checkbox" id="{id_}">'
                f'<div class="dot" style="background:{color};"></div>'
                f'<span>{label} ({count})</span></label>\n')

    chk_disease = _chk("chk_disease", "#e63946", "Disease variants", len(disease_vars))
    chk_other   = _chk("chk_other",   "#f4a261", "Other variants",   len(other_vars))
    chk_disulf  = _chk("chk_disulf",  "#FFD700", "Disulfide bonds",  len(disulfide_bonds))
    chk_glycan  = _chk("chk_glycan",  "#2dc653", "Glycosylation",    len(glycosylation))
    chk_active  = _chk("chk_active",  "#c77dff", "Active sites",     len(active_sites))
    chk_binding = _chk("chk_binding", "#4cc9f0", "Binding sites",    len(binding_sites))

    def _leg(color: str, label: str) -> str:
        return (f'<div class="leg-row"><div class="leg-swatch" style="background:{color};"></div>'
                f' {label}</div>\n')

    extra_legend = ""
    if disulfide_bonds: extra_legend += _leg("#FFD700", "Disulfide bond")
    if glycosylation:   extra_legend += _leg("#2dc653", "Glycosylation site")
    if active_sites:    extra_legend += _leg("#c77dff", "Active site")
    if binding_sites:   extra_legend += _leg("#4cc9f0", "Binding site")

    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.4/jquery.min.js"></script>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  body {{ margin:0; padding:0; background:#0d1117; font-family:sans-serif; color:#cdd3de; }}
  #controls {{
    position:absolute; top:10px; left:10px; z-index:100;
    background:rgba(13,17,23,0.88); border-radius:8px; padding:10px 14px;
    font-size:13px; min-width:230px;
  }}
  #controls h4 {{ margin:0; color:#e6edf3; font-size:14px; }}
  .ctrl-header {{ display:flex; justify-content:space-between; align-items:center;
                  margin-bottom:8px; }}
  .ctrl-collapse {{ background:none; border:none; color:#8b949e; cursor:pointer;
                    font-size:13px; padding:0 2px; line-height:1; }}
  .ctrl-collapse:hover {{ color:#e6edf3; }}
  .ctrl-row {{ margin:3px 0; display:flex; align-items:center; gap:6px; cursor:pointer; }}
  .ctrl-section {{ color:#8b949e; font-size:10px; text-transform:uppercase;
                   letter-spacing:0.8px; margin:8px 0 4px; border-top:1px solid #21262d;
                   padding-top:6px; }}
  .ctrl-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:2px 6px; margin-bottom:2px; }}
  .dot {{ width:12px; height:12px; border-radius:50%; flex-shrink:0; }}
  .bg-row {{ display:flex; align-items:center; gap:5px; margin-top:5px; }}
  .bg-btn {{ width:18px; height:18px; border-radius:3px; border:2px solid transparent;
             cursor:pointer; flex-shrink:0; }}
  .bg-btn.active {{ border-color:#e6edf3; }}
  .ctrl-btn {{ background:rgba(255,255,255,0.08); border:1px solid #30363d; color:#cdd3de;
               border-radius:4px; padding:3px 9px; cursor:pointer; font-size:12px;
               transition:background 0.15s; }}
  .ctrl-btn:hover {{ background:rgba(255,255,255,0.18); }}
  #shot-toast {{ display:none; position:absolute; bottom:14px; right:14px; z-index:200;
                 background:rgba(46,160,67,0.9); color:white; border-radius:6px;
                 padding:5px 12px; font-size:12px; pointer-events:none; }}
  #legend {{
    position:absolute; bottom:10px; left:10px; z-index:100;
    background:rgba(13,17,23,0.88); border-radius:8px; padding:10px 14px; font-size:12px;
  }}
  #legend h4 {{ margin:0 0 6px; color:#e6edf3; font-size:13px; }}
  .leg-row {{ display:flex; align-items:center; gap:6px; margin:3px 0; }}
  .leg-swatch {{ width:14px; height:14px; border-radius:2px; flex-shrink:0; }}
  #viewer {{ width:{width}px; height:{height}px; position:relative; }}
  #seq-panel {{
    width:{width}px; box-sizing:border-box;
    background:#161b22; border:1px solid #21262d; border-bottom:none;
    padding:7px 12px; white-space:pre-wrap; word-break:break-all;
    font-family:monospace; font-size:13px; line-height:1.6; user-select:none;
  }}
  .seq-aa {{ cursor:pointer; padding:0 1px; border-radius:2px; transition:background 0.1s; }}
  .seq-aa:hover {{ background:#30363d; }}
  .seq-sep {{ color:#8b949e; cursor:default; }}
  .seq-chain-lbl {{ color:#8b949e; font-size:11px; font-family:sans-serif; font-weight:600;
                    cursor:default; letter-spacing:0.5px; }}
</style>
</head><body>
<div id="seq-panel"><span id="seq-letters"></span></div>
<div style="position:relative;width:{width}px;height:{height}px;">
  <div id="viewer"></div>

  <div id="controls">
    <div class="ctrl-header">
      <h4>🧬 {protein_name}</h4>
      <button class="ctrl-collapse" onclick="toggleControls()" title="Collapse/expand" id="ctrl-toggle">▲</button>
    </div>
    <div id="ctrl-body">
    <div class="ctrl-section">Style</div>
    <div class="ctrl-grid">
      <label class="ctrl-row"><input type="radio" name="style" value="cartoon" checked><span>Cartoon</span></label>
      <label class="ctrl-row"><input type="radio" name="style" value="stick"><span>Stick</span></label>
      <label class="ctrl-row"><input type="radio" name="style" value="sphere"><span>Sphere</span></label>
      <label class="ctrl-row"><input type="radio" name="style" value="line"><span>Line</span></label>
      <label class="ctrl-row"><input type="radio" name="style" value="surface"><span>Surface</span></label>
    </div>
    <div id="surf-opacity-row" style="display:none;margin:4px 0 2px;align-items:center;gap:6px;font-size:12px;color:#8b949e;">
      <span>Opacity</span>
      <input type="range" id="surf-opacity" min="0" max="1" step="0.05" value="0.85"
             style="flex:1;accent-color:#58a6ff;" oninput="updateStyle()">
      <span id="surf-opacity-val">0.85</span>
    </div>

    <div class="ctrl-section">Color</div>
    <div class="ctrl-grid">
      <label class="ctrl-row"><input type="radio" name="color" value="plddt" checked><span>pLDDT</span></label>
      <label class="ctrl-row"><input type="radio" name="color" value="spectrum"><span>Spectrum</span></label>
      <label class="ctrl-row"><input type="radio" name="color" value="ss"><span>Sec. struct.</span></label>
      <label class="ctrl-row"><input type="radio" name="color" value="chain"><span>Chain</span></label>
    </div>

    <div class="ctrl-section">Annotations</div>
    {chk_disease}{chk_other}{chk_disulf}{chk_glycan}{chk_active}{chk_binding}
    <label class="ctrl-row"><input type="checkbox" id="chk_labels"><span>Variant labels</span></label>

    <div class="ctrl-section">Options</div>
    <label class="ctrl-row"><input type="checkbox" id="chk_spin"><span>Spin</span></label>
    <label class="ctrl-row"><input type="checkbox" id="chk_ortho"><span>Orthographic</span></label>
    <div class="bg-row">
      <span style="color:#8b949e;font-size:11px;">BG</span>
      <div class="bg-btn active" id="bg0" style="background:#0d1117;" title="Dark" onclick="setBg('#0d1117','bg0')"></div>
      <div class="bg-btn" id="bg1" style="background:#000000;" title="Black" onclick="setBg('#000000','bg1')"></div>
      <div class="bg-btn" id="bg2" style="background:#ffffff;" title="White" onclick="setBg('#ffffff','bg2')"></div>
      <div class="bg-btn" id="bg3" style="background:#3d444d;" title="Grey" onclick="setBg('#3d444d','bg3')"></div>
    </div>

    <div class="ctrl-section">Screenshot</div>
    <label class="ctrl-row"><input type="checkbox" id="chk_transparent"><span>Transparent BG</span></label>
    <div style="display:flex;gap:5px;margin-top:5px;">
      <button class="ctrl-btn" onclick="screenshotDownload()" title="Download PNG">📥 Save</button>
      <button class="ctrl-btn" onclick="screenshotClipboard()" title="Copy to clipboard">📋 Copy</button>
    </div>
    </div><!-- /ctrl-body -->
  </div>

  <div id="legend">
    <div id="legend-plddt">
      <h4>pLDDT confidence</h4>
      <div class="leg-row"><div class="leg-swatch" style="background:#0053d6;"></div> ≥90 Very high</div>
      <div class="leg-row"><div class="leg-swatch" style="background:#65cbf3;"></div> 70–90 High</div>
      <div class="leg-row"><div class="leg-swatch" style="background:#ffdb13;"></div> 50–70 Low</div>
      <div class="leg-row"><div class="leg-swatch" style="background:#ff7d45;"></div> &lt;50 Very low</div>
    </div>
    <div id="legend-other" style="display:none;"><h4 id="legend-other-title"></h4></div>
    {extra_legend}
  </div>

  <div id="variant-panel" style="
    display:none; position:absolute; top:10px; right:10px; z-index:100;
    background:rgba(13,17,23,0.93); border:1px solid rgba(230,57,70,0.55);
    border-radius:8px; padding:12px 14px; width:250px;
    font-size:13px; box-shadow:0 4px 24px rgba(0,0,0,0.55);">
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;">
      <span id="vp-mutation" style="font-size:22px;font-weight:700;color:#ff6b81;letter-spacing:0.5px;"></span>
      <button onclick="closeVariantPanel()" style="background:rgba(255,255,255,0.08);border:none;
        color:#cdd3de;cursor:pointer;border-radius:4px;width:22px;height:22px;
        font-size:14px;line-height:1;flex-shrink:0;">✕</button>
    </div>
    <div style="font-size:10px;text-transform:uppercase;letter-spacing:0.7px;color:#8b949e;margin-bottom:3px;">Disease / Condition</div>
    <div id="vp-disease" style="color:#e6edf3;margin-bottom:9px;font-weight:500;"></div>
    <div style="font-size:10px;text-transform:uppercase;letter-spacing:0.7px;color:#8b949e;margin-bottom:3px;">Description</div>
    <div id="vp-desc" style="color:#cdd3de;font-size:12px;line-height:1.45;margin-bottom:7px;"></div>
    <div id="vp-pos" style="color:#8b949e;font-size:11px;"></div>
  </div>
  <div id="shot-toast">✓ Copied to clipboard</div>
</div>

<script>
var pdbData        = `{pdb_str}`;
var diseaseVars    = {disease_js};
var otherVars      = {other_js};
var disulfideBonds = {disulfide_js};
var glycanSites    = {glycan_js};
var activeSites    = {active_js};
var bindingSites   = {binding_js};
var seqData        = {seq_js};

$(document).ready(function() {{
  var viewer = $3Dmol.createViewer('viewer', {{backgroundColor:'#0d1117', antialias:true}});
  viewer.addModel(pdbData, 'pdb');

  var selectedResi  = null;
  var selectedChain = null;
  var labelHandles  = [];
  var disulfLines   = [];
  var currentSurface = null;
  var surfaceEpoch   = 0;   // incremented each renderAll; lets late callbacks self-cancel

  function getColorSpec() {{
    var s = $('input[name="color"]:checked').val();
    if (s === 'plddt') return {{colorfunc: function(atom) {{
      var b=atom.b;
      if(b>=90) return '#0053d6';
      if(b>=70) return '#65cbf3';
      if(b>=50) return '#ffdb13';
      return '#ff7d45';
    }}}};
    if (s === 'spectrum') return {{colorscheme:'spectrum'}};
    if (s === 'ss')       return {{colorscheme:'ssJmol'}};
    if (s === 'chain')    return {{colorscheme:'chainHetatm'}};
    return {{colorscheme:'Jmol'}};
  }}

  function updateLegend() {{
    var s = $('input[name="color"]:checked').val();
    if (s === 'plddt') {{
      $('#legend-plddt').show(); $('#legend-other').hide();
    }} else {{
      $('#legend-plddt').hide();
      var titles = {{spectrum:'N→C spectrum',ss:'Secondary structure',chain:'By chain'}};
      $('#legend-other-title').text(titles[s]||s);
      $('#legend-other').show();
    }}
  }}

  var currentBg = '#0d1117';

  window.toggleControls = function() {{
    var body = document.getElementById('ctrl-body');
    var btn = document.getElementById('ctrl-toggle');
    var collapsed = body.style.display === 'none';
    body.style.display = collapsed ? '' : 'none';
    btn.textContent = collapsed ? '▲' : '▼';
  }};

  window.setBg = function(color, btnId) {{
    currentBg = color;
    viewer.setBackgroundColor(color); viewer.render();
    $('.bg-btn').removeClass('active');
    $('#'+btnId).addClass('active');
  }};

  function parseBgRGB(hex) {{
    var c = document.createElement('canvas'); c.width = c.height = 1;
    var ctx = c.getContext('2d');
    ctx.fillStyle = hex; ctx.fillRect(0,0,1,1);
    return ctx.getImageData(0,0,1,1).data;
  }}

  function captureURI() {{
    var baseURI = viewer.pngURI();
    if (!$('#chk_transparent').prop('checked')) return Promise.resolve(baseURI);
    return new Promise(function(resolve) {{
      var img = new Image();
      img.onload = function() {{
        var c = document.createElement('canvas');
        c.width = img.width; c.height = img.height;
        var ctx = c.getContext('2d');
        ctx.drawImage(img, 0, 0);
        var d = ctx.getImageData(0, 0, c.width, c.height);
        var px = d.data;
        var bg = parseBgRGB(currentBg);
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

  window.screenshotDownload = function() {{
    captureURI().then(function(uri) {{
      var a = document.createElement('a');
      a.href = uri; a.download = 'structure.png'; a.click();
    }});
  }};

  window.screenshotClipboard = function() {{
    captureURI().then(function(uri) {{
      fetch(uri).then(function(r) {{ return r.blob(); }}).then(function(blob) {{
        navigator.clipboard.write([new ClipboardItem({{'image/png': blob}})])
          .then(function() {{
            var t = document.getElementById('shot-toast');
            t.style.display = 'block';
            setTimeout(function() {{ t.style.display = 'none'; }}, 1800);
          }})
          .catch(function() {{
            alert('Clipboard blocked by browser — use Save, or open in Chrome/Edge over HTTPS.');
          }});
      }});
    }});
  }};

  // Precompute feature lookup sets once
  var dSet={{}},oSet={{}},aSet={{}},bSet={{}},gSet={{}},sSet={{}};
  diseaseVars.forEach(function(v)    {{ dSet[v.resi]=1; }});
  otherVars.forEach(function(v)      {{ oSet[v.resi]=1; }});
  activeSites.forEach(function(a)    {{ aSet[a.resi]=1; }});
  bindingSites.forEach(function(b)   {{ bSet[b.resi]=1; }});
  glycanSites.forEach(function(g)    {{ gSet[g.resi]=1; }});
  disulfideBonds.forEach(function(b) {{ sSet[b.pos1]=1; sSet[b.pos2]=1; }});

  // Build sequence strip + cache DOM refs
  var multiChain = seqData.length > 1 && seqData.some(function(r) {{ return r.chain !== seqData[0].chain; }});
  var seqEls = {{}};
  (function() {{
    var html = '';
    var prevChain = null;
    var chainPos  = 0;
    seqData.forEach(function(r) {{
      if (r.chain !== prevChain) {{
        if (prevChain !== null) html += '<br>';
        if (multiChain) html += '<span class="seq-chain-lbl">Chain ' + r.chain + '</span><br>';
        prevChain = r.chain;
        chainPos  = 0;
      }}
      html += '<span class="seq-aa" id="seq-'+r.chain+'-'+r.resi+'" data-resi="'+r.resi
            + '" data-chain="'+r.chain+'" title="'+r.resn+' '+r.resi+'">'+r.aa+'</span>';
      chainPos++;
      if      (chainPos % 100 === 0) html += '<br>';
      else if (chainPos % 50  === 0) html += '<span class="seq-sep">  </span>';
      else if (chainPos % 10  === 0) html += '<span class="seq-sep"> </span>';
    }});
    document.getElementById('seq-letters').innerHTML = html;
    seqData.forEach(function(r) {{ seqEls[r.chain+'-'+r.resi]=document.getElementById('seq-'+r.chain+'-'+r.resi); }});
  }})();

  function renderSeq() {{
    var showD=document.getElementById('chk_disease')&&document.getElementById('chk_disease').checked;
    var showO=document.getElementById('chk_other')  &&document.getElementById('chk_other').checked;
    var showA=document.getElementById('chk_active') &&document.getElementById('chk_active').checked;
    var showB=document.getElementById('chk_binding')&&document.getElementById('chk_binding').checked;
    var showG=document.getElementById('chk_glycan') &&document.getElementById('chk_glycan').checked;
    var showS=document.getElementById('chk_disulf') &&document.getElementById('chk_disulf').checked;
    seqData.forEach(function(r) {{
      var el=seqEls[r.chain+'-'+r.resi]; if(!el) return;
      if(r.resi===selectedResi&&r.chain===selectedChain) {{
        el.style.background='#ff1744'; el.style.color='white'; el.style.borderRadius='2px'; return;
      }}
      var c='#566070';
      if      (showD&&dSet[r.resi]) {{ c='#e63946'; }}
      else if (showO&&oSet[r.resi]) {{ c='#f4a261'; }}
      else if (showA&&aSet[r.resi]) {{ c='#c77dff'; }}
      else if (showB&&bSet[r.resi]) {{ c='#4cc9f0'; }}
      else if (showG&&gSet[r.resi]) {{ c='#2dc653'; }}
      else if (showS&&sSet[r.resi]) {{ c='#FFD700'; }}
      el.style.background=''; el.style.color=c; el.style.borderRadius='';
    }});
    var selKey = selectedChain+'-'+selectedResi;
    if(selectedResi!==null&&seqEls[selKey])
      seqEls[selKey].scrollIntoView({{behavior:'smooth',block:'nearest',inline:'center'}});
  }}

  function renderAll() {{
    if (currentSurface !== null) {{
      viewer.removeSurface(currentSurface); currentSurface = null;
    }}
    surfaceEpoch++;
    var myEpoch = surfaceEpoch;

    viewer.setStyle({{}},{{}});

    var styleVal = $('input[name="style"]:checked').val();
    var col      = getColorSpec();
    var opacityRow = document.getElementById('surf-opacity-row');
    opacityRow.style.display = (styleVal === 'surface') ? 'flex' : 'none';
    if      (styleVal === 'cartoon') viewer.setStyle({{}},{{cartoon:col}});
    else if (styleVal === 'stick')   viewer.setStyle({{}},{{stick:Object.assign({{radius:0.15}},col)}});
    else if (styleVal === 'sphere')  viewer.setStyle({{}},{{sphere:Object.assign({{scale:0.4}},col)}});
    else if (styleVal === 'line')    viewer.setStyle({{}},{{line:col}});
    else if (styleVal === 'surface') {{
      var surfOpacity = parseFloat(document.getElementById('surf-opacity').value);
      document.getElementById('surf-opacity-val').textContent = surfOpacity.toFixed(2);
      viewer.setStyle({{}},{{cartoon:{{opacity:0.25,color:'#888888'}}}});
      var surfP = viewer.addSurface($3Dmol.SurfaceType.VDW, Object.assign({{opacity:surfOpacity}},col), {{hetflag:false}});
      Promise.resolve(surfP).then(function(id) {{
        if (myEpoch !== surfaceEpoch) {{
          viewer.removeSurface(id);  // style changed while computing — discard
        }} else {{
          currentSurface = id;
        }}
      }});
    }}

    if($('#chk_disease').prop('checked'))
      diseaseVars.forEach(function(v)  {{ viewer.addStyle({{resi:v.resi}},{{stick:{{radius:0.2,color:'#e63946'}}}}); }});
    if($('#chk_other').prop('checked'))
      otherVars.forEach(function(v)    {{ viewer.addStyle({{resi:v.resi}},{{stick:{{radius:0.18,color:'#f4a261'}}}}); }});
    if($('#chk_disulf').prop('checked'))
      disulfideBonds.forEach(function(b) {{
        viewer.addStyle({{resi:b.pos1}},{{stick:{{radius:0.2,color:'#FFD700'}}}});
        viewer.addStyle({{resi:b.pos2}},{{stick:{{radius:0.2,color:'#FFD700'}}}});
      }});
    if($('#chk_glycan').prop('checked'))
      glycanSites.forEach(function(g)  {{ viewer.addStyle({{resi:g.resi}},{{stick:{{radius:0.18,color:'#2dc653'}}}}); }});
    if($('#chk_active').prop('checked'))
      activeSites.forEach(function(a)  {{ viewer.addStyle({{resi:a.resi}},{{stick:{{radius:0.18,color:'#c77dff'}}}}); }});
    if($('#chk_binding').prop('checked'))
      bindingSites.forEach(function(b) {{ viewer.addStyle({{resi:b.resi}},{{stick:{{radius:0.18,color:'#4cc9f0'}}}}); }});

    if(selectedResi!==null) {{
      viewer.addStyle({{resi:selectedResi}},{{stick:{{radius:0.26,colorscheme:'Jmol'}}}});
      viewer.addStyle({{resi:selectedResi,atom:'CA'}},{{sphere:{{color:'#ff1744',radius:0.45}}}});
    }}

    diseaseVars.forEach(function(v) {{
      viewer.setClickable({{resi:v.resi}},true,function(){{ selectVariant(v); }});
    }});

    updateLegend();
    renderSeq();
    viewer.render();
  }}

  function drawDisulfLines() {{
    disulfLines.forEach(function(ln) {{ viewer.removeLine(ln); }});
    disulfLines=[];
    if(!$('#chk_disulf').prop('checked')) return;
    disulfideBonds.forEach(function(b) {{
      var sg1=viewer.getModel().selectedAtoms({{resi:b.pos1,atom:'SG'}});
      var sg2=viewer.getModel().selectedAtoms({{resi:b.pos2,atom:'SG'}});
      if(sg1.length&&sg2.length)
        disulfLines.push(viewer.addLine({{
          start:{{x:sg1[0].x,y:sg1[0].y,z:sg1[0].z}},
          end:  {{x:sg2[0].x,y:sg2[0].y,z:sg2[0].z}},
          linewidth:3,color:'#FFD700'
        }}));
    }});
  }}

  function drawLabels() {{
    labelHandles.forEach(function(h) {{ viewer.removeLabel(h); }});
    labelHandles=[];
    if(!$('#chk_labels').prop('checked')) return;
    diseaseVars.forEach(function(v) {{
      labelHandles.push(viewer.addLabel(v.label,{{
        resi:v.resi,atom:'CA',fontSize:11,fontColor:'white',
        backgroundColor:'rgba(230,57,70,0.75)',
        borderColor:'#e63946',borderThickness:1,padding:2,inFront:true
      }}));
    }});
  }}

  function showVariantPanel(v) {{
    $('#vp-mutation').text(v.label);
    $('#vp-disease').text(v.disease||'—');
    $('#vp-desc').text(v.desc||'No description available.');
    $('#vp-pos').text('Residue position: '+v.resi);
    $('#variant-panel').show();
  }}

  window.closeVariantPanel = function() {{
    selectedResi=null; selectedChain=null;
    $('#variant-panel').hide();
    viewer.zoomTo();
    renderAll();
  }};

  function selectVariant(v) {{
    selectedResi=v.resi; selectedChain=seqData.length>0?seqData[0].chain:null;
    viewer.center({{resi:v.resi}},800,false);
    renderAll();
    showVariantPanel(v);
  }}

  $(document).on('click','.seq-aa',function() {{
    var resi=parseInt($(this).data('resi'));
    var chain=$(this).data('chain');
    selectedResi=resi; selectedChain=chain;
    viewer.center({{resi:resi}},500,false);
    var dv=diseaseVars.find(function(v){{ return v.resi===resi; }});
    if(dv) {{ showVariantPanel(dv); }} else {{ $('#variant-panel').hide(); }}
    renderAll();
  }});

  renderAll();
  drawDisulfLines();
  drawLabels();
  viewer.zoomTo();
  viewer.render();

  $('input[name="style"]').change(function() {{ renderAll(); }});
  $('input[name="color"]').change(function() {{ renderAll(); }});
  $('#chk_disease').change(function(){{ renderAll(); }});
  $('#chk_other').change(function()  {{ renderAll(); }});
  $('#chk_glycan').change(function() {{ renderAll(); }});
  $('#chk_active').change(function() {{ renderAll(); }});
  $('#chk_binding').change(function(){{ renderAll(); }});
  $('#chk_disulf').change(function() {{ drawDisulfLines(); renderAll(); }});
  $('#chk_labels').change(function() {{ drawLabels(); viewer.render(); }});
  $('#chk_spin').change(function() {{
    $(this).prop('checked') ? viewer.spin('y',1) : viewer.spin(false);
  }});
  $('#chk_ortho').change(function() {{
    viewer.setProjection($(this).prop('checked') ? 'orthographic' : 'perspective');
    viewer.render();
  }});
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
        help="Standard PDB format files (.pdb or .ent)",
        label_visibility="collapsed",
    )
with opt_col:
    with st.popover("⚙", help="Analysis settings"):
        contact_cutoff = st.slider(
            "Contact cutoff (Å)", min_value=3.0, max_value=7.0,
            value=4.5, step=0.1,
            help="Maximum heavy-atom distance for non-bonded contacts",
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

    st.title(f"🧬 {pdb_path.stem.upper()}")
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
    c2.metric("R-work",     str(summary["r_work"]) if summary["r_work"] else "N/A")
    c3.metric("R-free",     str(summary["r_free"]) if summary["r_free"] else "N/A")
    c4.metric("Residues",   f"{summary['n_residues']:,}")
    c5.metric("Ligands",    summary["n_ligands"])
    c6.metric("Rg (Å)",     f"{rg:.2f}" if rg else "N/A")
    st.divider()

else:
    pdb_bytes = pdb_path = summary = bf_data = rama_data = None
    rg = fasta_str = missing_data = reports = None


# ═══════════════════════════════════════════════════════════════════════════════
# Tabs — always rendered; AlphaFold tab works without a PDB
# ═══════════════════════════════════════════════════════════════════════════════
tab_summary, tab_bsa, tab_viewer, tab_2d, tab_af = st.tabs([
    "📋 Structure Summary",
    "⬛ Buried Surface Area",
    "🌐 3D Viewer",
    "🖼 2D Diagrams",
    "🧬 AlphaFold + Variants",
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
        bsa_head_col, bsa_opt_col = st.columns([11, 1])
        with bsa_opt_col:
            with st.popover("⚙", help="BSA settings"):
                probe_radius = st.slider(
                    "Probe radius (Å)", min_value=1.0, max_value=2.0,
                    value=1.4, step=0.1,
                    help="Rolling sphere radius; 1.4 Å = water molecule",
                )
                n_points_bsa = st.select_slider(
                    "Sphere points (accuracy)", options=[50, 100, 200, 500],
                    value=100,
                    help="More points = higher accuracy, slower",
                )
        st.caption("BSA(A, B) = ( SASA_A + SASA_B − SASA_AB ) / 2")
        n_chains = len(summary["chains"].split(",")) if summary["chains"] else 0
        if n_chains < 2:
            st.warning("Only one chain detected — BSA requires at least two chains.")
        else:
            with st.spinner(f"Computing BSA for all chain pairs (probe={probe_radius} Å)…"):
                bsa_results, chain_ids = _run_bsa(pdb_bytes, probe_radius, n_points_bsa)
            if not bsa_results:
                st.info("No inter-chain contacts found.")
            else:
                top = bsa_results[0]
                n_iface = sum(1 for r in bsa_results if r["bsa_total"] > 0)
                st.success(
                    f"**{n_iface}** interface(s) detected across **{len(chain_ids)}** chains  ·  "
                    f"Largest: chains {top['chain_1']}–{top['chain_2']} ({top['bsa_total']:.0f} Å²)"
                )
                bsa_png = _png_bsa_matrix(bsa_results, chain_ids, title=pdb_path.stem.upper())
                st.image(bsa_png, width="stretch")
                st.download_button("⬇ BSA matrix PNG", bsa_png,
                                   file_name=f"{pdb_path.stem}_bsa_matrix.png", mime="image/png")
                st.divider()
                st.markdown("#### Pairwise BSA table")
                df_bsa = pd.DataFrame(bsa_results)
                df_bsa.columns = [
                    "Chain 1", "Chain 2",
                    "SASA 1 (Å²)", "SASA 2 (Å²)", "SASA complex (Å²)",
                    "BSA total (Å²)", "BSA on 1 (Å²)", "BSA on 2 (Å²)",
                    "Interface % 1", "Interface % 2",
                ]
                st.dataframe(
                    df_bsa.style.background_gradient(subset=["BSA total (Å²)"], cmap="YlOrRd"),
                    use_container_width=True, hide_index=True,
                )
                st.download_button("⬇ BSA table CSV", df_bsa.to_csv(index=False),
                                   file_name=f"{pdb_path.stem}_bsa.csv", mime="text/csv")
                st.divider()
                st.markdown("#### Per-chain SASA summary")
                seen: dict = {}
                for r in bsa_results:
                    if r["chain_1"] not in seen: seen[r["chain_1"]] = r["sasa_1"]
                    if r["chain_2"] not in seen: seen[r["chain_2"]] = r["sasa_2"]
                for cid in chain_ids:
                    if cid not in seen: seen[cid] = None
                sasa_rows = []
                for cid in chain_ids:
                    sv = seen.get(cid)
                    tb = sum(
                        r["bsa_on_1"] if r["chain_1"] == cid else r["bsa_on_2"]
                        for r in bsa_results
                        if (r["chain_1"] == cid or r["chain_2"] == cid) and r["bsa_total"] > 0
                    )
                    sasa_rows.append({
                        "Chain": cid, "SASA isolated (Å²)": sv,
                        "Total buried (Å²)": round(tb, 1),
                        "% buried": round(100 * tb / sv, 1) if sv else None,
                    })
                st.dataframe(pd.DataFrame(sasa_rows), use_container_width=True, hide_index=True)


# ─── TAB 3: 3D Viewer + Ligand Analysis ───────────────────────────────────────
with tab_viewer:
    if not uploaded:
        st.info(_NO_PDB)
    else:
        st.subheader("Interactive 3D Viewer")
        if not reports:
            st.warning("No primary ligands detected — 3D viewer unavailable.")
        else:
            lig_labels = [f"{r.ligand_resn}  chain {r.ligand_chain}:{r.ligand_resi}" for r in reports]
            sel_idx = st.selectbox("Select ligand to visualise", range(len(reports)),
                                   format_func=lambda i: lig_labels[i]) if len(reports) > 1 else 0
            if len(reports) == 1:
                st.caption(f"Showing: **{lig_labels[0]}**")
            report = reports[sel_idx]
            st.caption(
                f"{report.summary()}  ·  "
                "Use toggles in the viewer panel to show/hide interaction types. "
                "Drag to rotate · scroll to zoom."
            )
            lc1, lc2, lc3, lc4, lc5 = st.columns(5)
            lc1.markdown("🟡 **Ligand**"); lc2.markdown("🔵 **Contacts**")
            lc3.markdown("🟠 **H-bonds**"); lc4.markdown("🟣 **Pi**")
            lc5.markdown("🔵🔴 **Salt bridges**")

            @st.cache_resource
            def _cached_structure(path_str: str):
                return al_load(path_str)

            structure = _cached_structure(str(pdb_path))
            html_str = _render_html(pdb_path, report, structure, width=1100, height=680)
            components.html(html_str, height=720, scrolling=False)
            st.download_button(
                "⬇ Download standalone HTML viewer", html_str,
                file_name=f"{report.ligand_resn}_{report.ligand_chain}{report.ligand_resi}.html",
                mime="text/html",
            )

        st.divider()
        st.subheader("Ligand Interaction Analysis")
        st.caption(f"Contact cutoff: {contact_cutoff} Å  ·  Adjust the slider at the top of the page to recompute")
        if not reports:
            st.warning(
                "No primary ligands detected.  \n"
                "All HETATM residues are water, crystallographic excipients, "
                "or have fewer than 7 heavy atoms."
            )
        else:
            st.success(f"**{len(reports)}** primary ligand(s) detected and analysed.")
            for rep in reports:
                lig_id = f"{rep.ligand_resn} chain {rep.ligand_chain} resi {rep.ligand_resi}"
                counts = (f"{rep.n_contacts} contacts · {rep.n_hbonds} H-bonds · "
                          f"{rep.n_pi} pi · {rep.n_salt_bridges} salt bridges")
                with st.expander(f"**{lig_id}**  —  {counts}", expanded=False):
                    it1, it2, it3, it4 = st.tabs([
                        f"Contacts ({rep.n_contacts})", f"H-bonds ({rep.n_hbonds})",
                        f"Pi interactions ({rep.n_pi})", f"Salt bridges ({rep.n_salt_bridges})",
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
                            st.download_button("⬇ CSV", df.to_csv(index=False),
                                               file_name=f"{rep.ligand_resn}_hbonds.csv",
                                               mime="text/csv", key=f"h_{key}")
                        else:
                            st.info("No hydrogen bonds detected.")
                    with it3:
                        if rep.pi_interactions:
                            df = pd.DataFrame([p.as_dict() for p in rep.pi_interactions])
                            st.dataframe(df, use_container_width=True, hide_index=True)
                            st.download_button("⬇ CSV", df.to_csv(index=False),
                                               file_name=f"{rep.ligand_resn}_pi.csv",
                                               mime="text/csv", key=f"p_{key}")
                        else:
                            st.info("No pi interactions detected.")
                    with it4:
                        if rep.salt_bridges:
                            df = pd.DataFrame([s.as_dict() for s in rep.salt_bridges])
                            st.dataframe(df, use_container_width=True, hide_index=True)
                            st.download_button("⬇ CSV", df.to_csv(index=False),
                                               file_name=f"{rep.ligand_resn}_salt.csv",
                                               mime="text/csv", key=f"s_{key}")
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
                st.caption(rep.summary())
                col_fp, col_2d = st.columns(2, gap="large")
                key = f"{rep.ligand_resn}_{rep.ligand_resi}"
                has_interactions = (rep.n_contacts + rep.n_hbonds +
                                    rep.n_pi + rep.n_salt_bridges) > 0
                with col_fp:
                    st.markdown("#### Interaction Fingerprint")
                    st.caption("Dot size = contact count · Dot shade = proximity · Number = count when >1")
                    if has_interactions:
                        fp_bytes = _png_fingerprint(rep)
                        st.image(fp_bytes, width="stretch")
                        st.download_button("⬇ Fingerprint PNG", fp_bytes,
                                           file_name=f"{rep.ligand_resn}_fingerprint.png",
                                           mime="image/png", key=f"fp_{key}")
                    else:
                        st.info("No interactions to plot.")
                with col_2d:
                    st.markdown("#### LIGPLOT-style 2D Diagram")
                    st.caption("Ligand in centre (CPK colours) · Protein residues at periphery · "
                               "Spoked arcs = hydrophobic contacts")
                    if has_interactions:
                        d2_bytes = _png_2d(pdb_path, rep, structure_2d)
                        st.image(d2_bytes, width="stretch")
                        st.download_button("⬇ 2D diagram PNG", d2_bytes,
                                           file_name=f"{rep.ligand_resn}_2d.png",
                                           mime="image/png", key=f"2d_{key}")
                    else:
                        st.info("No interactions to plot.")
                st.divider()


# ─── TAB 5: AlphaFold + Variants ──────────────────────────────────────────────
with tab_af:
    st.subheader("AlphaFold Structure + Disease Variants")
    st.caption(
        "Fetches the AlphaFold predicted structure and disease-associated SNPs "
        "from UniProt. No PDB file required."
    )

    af_uid = st.text_input(
        "UniProt accession",
        placeholder="e.g. P68871 (HBB), P04637 (TP53)",
        help="Enter a UniProt accession ID to fetch the AlphaFold model and variants.",
    ).strip()

    if af_uid:
        with st.spinner(f"Fetching AlphaFold structure for {af_uid}…"):
            af_pdb_bytes = _fetch_alphafold_pdb(af_uid)

        if af_pdb_bytes is None:
            st.error(
                f"Could not fetch AlphaFold structure for **{af_uid}**. "
                "Check the accession and try again."
            )
        else:
            with st.spinner(f"Loading UniProt annotations for {af_uid}…"):
                uniprot_data  = _fetch_uniprot(af_uid)
                humsavar_db   = _load_humsavar()
                all_variants  = humsavar_db.get(af_uid.upper(), [])

            features: dict = {}
            gene = af_uid
            if uniprot_data:
                features = _parse_uniprot_features(uniprot_data)
                genes = uniprot_data.get("genes", [])
                if genes:
                    gene = genes[0].get("geneName", {}).get("value", af_uid)

            disease_vars = [v for v in all_variants if v["is_disease"]]
            other_vars   = [v for v in all_variants if not v["is_disease"]]

            mc1, mc2, mc3, mc4 = st.columns(4)
            mc1.metric("Disease variants", len(disease_vars))
            mc2.metric("Other variants",   len(other_vars))
            mc3.metric("Disulfide bonds",  len(features.get("disulfide_bonds", [])))
            mc4.metric("Active/Binding sites",
                       len(features.get("active_sites", [])) +
                       len(features.get("binding_sites", [])))

            feat_counts = (
                f"{len(features.get('glycosylation', []))} glycosylation  ·  "
                f"{len(features.get('active_sites', []))} active sites  ·  "
                f"{len(features.get('binding_sites', []))} binding sites"
            )
            st.caption(
                f"Coloured by pLDDT confidence · {feat_counts}  ·  "
                "Toggle layers in the viewer panel"
            )

            af_html = _render_alphafold_html(
                af_pdb_bytes, all_variants, features=features,
                protein_name=f"{gene} — {af_uid}",
                width=1100, height=680,
            )
            components.html(af_html, height=820, scrolling=False)

            dl_c1, dl_c2 = st.columns(2)
            with dl_c1:
                st.download_button(
                    "⬇ Download AlphaFold PDB", af_pdb_bytes,
                    file_name=f"AF-{af_uid}-F1-model.pdb",
                    mime="chemical/x-pdb",
                )

            st.divider()

            if disease_vars:
                st.markdown("#### Disease-associated variants")
                df_dv = pd.DataFrame([{
                    "Position":    v["position"],
                    "Mutation":    v["mutation"],
                    "Disease":     v["disease"],
                    "Description": v["description"],
                } for v in disease_vars])
                st.dataframe(df_dv, use_container_width=True, hide_index=True)
                st.download_button(
                    "⬇ Disease variants CSV", df_dv.to_csv(index=False),
                    file_name=f"{af_uid}_disease_variants.csv", mime="text/csv",
                )

            if other_vars:
                with st.expander(f"Other variants ({len(other_vars)})"):
                    df_ov = pd.DataFrame([{
                        "Position":    v["position"],
                        "Mutation":    v["mutation"],
                        "Description": v["description"],
                    } for v in other_vars])
                    st.dataframe(df_ov, use_container_width=True, hide_index=True)
                    st.download_button(
                        "⬇ Other variants CSV", df_ov.to_csv(index=False),
                        file_name=f"{af_uid}_other_variants.csv", mime="text/csv",
                    )
    else:
        st.info(
            "Enter a UniProt accession above to load the AlphaFold predicted structure "
            "and annotate disease variants, disulfide bonds, glycosylation, "
            "active sites, and binding sites."
        )
