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
    initial_sidebar_state="expanded",
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

_DISEASE_KEYWORDS = {
    "disease", "syndrome", "deficiency", "disorder", "cancer", "carcinoma",
    "tumor", "tumour", "dystrophy", "dysplasia", "neuropathy", "myopathy",
    "ataxia", "fibrosis", "retardation", "epilepsy", "diabetic", "diabetes",
    "hypertension", "cardiomyopathy", "anemia", "anaemia", "thrombosis",
    "leukemia", "leukaemia", "lymphoma", "melanoma", "mutation associated",
    "pathogenic", "thalassemia", "thalassaemia", "sickle", "polycythemia",
    "polycythaemia", "methemoglobin", "methaemoglobin", "hemolytic", "haemolytic",
}
_OMIM_CODE_RE = re.compile(
    r"^[Ii]n\s+[A-Z][A-Z0-9]*(?:-[A-Z0-9]+)+\s*[;,]"
    r"|^[Ii]n\s+[A-Z]{2,}[A-Z0-9-]*\s*[;,]"
)


@st.cache_data(show_spinner=False, ttl=3600)
def _fetch_alphafold_pdb(uniprot_id: str) -> bytes | None:
    uid = uniprot_id.upper()
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
            return r.read()
    except Exception:
        return None


@st.cache_data(show_spinner=False, ttl=3600)
def _fetch_uniprot(uniprot_id: str) -> dict | None:
    uid = uniprot_id.upper()
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=20) as r:
            return json.loads(r.read())
    except Exception:
        return None


def _parse_variants(data: dict) -> list[dict]:
    variants = []
    for feat in data.get("features", []):
        if feat.get("type") != "Natural variant":
            continue
        loc   = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        if start is None:
            continue
        desc     = feat.get("description", "")
        orig     = feat.get("alternativeSequence", {}).get("originalSequence", "")
        alt_list = feat.get("alternativeSequence", {}).get("alternativeSequences", [])
        alt      = alt_list[0] if alt_list else ""
        mutation = f"{orig}{start}{alt}" if (orig and alt) else f"pos {start}"
        desc_lower = desc.lower()
        is_disease = bool(_OMIM_CODE_RE.match(desc)) or any(
            kw in desc_lower for kw in _DISEASE_KEYWORDS
        )
        disease_name = ""
        if is_disease:
            m = re.search(r"[Ii]n ([^;,]+)", desc)
            disease_name = m.group(1).strip() if m else desc[:60]
        xrefs = {x.get("database"): x.get("id") for x in feat.get("evidences", [])}
        variants.append({
            "position":    start,
            "mutation":    mutation,
            "description": desc,
            "disease":     disease_name,
            "is_disease":  is_disease,
            "dbsnp":       xrefs.get("dbSNP", ""),
            "clinvar":     xrefs.get("ClinVar", ""),
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

    # Extract sequence from CA atoms
    seq_residues: list[dict] = []
    seen_resi: set[int] = set()
    for line in pdb_bytes.decode("utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") and len(line) > 25 and line[12:16].strip() == "CA":
            try:
                resi = int(line[22:26].strip())
                resn = line[17:20].strip()
                if resi not in seen_resi:
                    seen_resi.add(resi)
                    seq_residues.append({"resi": resi, "resn": resn, "aa": _AA3TO1.get(resn, "?")})
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
    font-size:13px; min-width:215px;
  }}
  #controls h4 {{ margin:0 0 8px; color:#e6edf3; font-size:14px; }}
  .ctrl-row {{ margin:4px 0; display:flex; align-items:center; gap:6px; cursor:pointer; }}
  .dot {{ width:12px; height:12px; border-radius:50%; flex-shrink:0; }}
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
    background:#161b22; border:1px solid #21262d; border-top:none;
    padding:7px 12px; overflow-x:auto; white-space:nowrap;
    font-family:monospace; font-size:13px; line-height:1.6; user-select:none;
  }}
  .seq-lbl {{ color:#8b949e; font-size:11px; font-family:sans-serif; margin-right:6px; }}
  .seq-aa {{ cursor:pointer; padding:0 1px; border-radius:2px; transition:background 0.1s; }}
  .seq-aa:hover {{ background:#30363d; }}
</style>
</head><body>
<div style="position:relative;width:{width}px;height:{height}px;">
  <div id="viewer"></div>

  <div id="controls">
    <h4>🧬 {protein_name}</h4>
    <label class="ctrl-row"><input type="checkbox" id="chk_cartoon" checked><span>Protein (cartoon)</span></label>
    {chk_disease}{chk_other}{chk_disulf}{chk_glycan}{chk_active}{chk_binding}
    <label class="ctrl-row"><input type="checkbox" id="chk_labels"><span>Variant labels</span></label>
  </div>

  <div id="legend">
    <h4>pLDDT confidence</h4>
    <div class="leg-row"><div class="leg-swatch" style="background:#0053d6;"></div> ≥90 Very high</div>
    <div class="leg-row"><div class="leg-swatch" style="background:#65cbf3;"></div> 70–90 High</div>
    <div class="leg-row"><div class="leg-swatch" style="background:#ffdb13;"></div> 50–70 Low</div>
    <div class="leg-row"><div class="leg-swatch" style="background:#ff7d45;"></div> &lt;50 Very low</div>
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
</div>

<div id="seq-panel"><span class="seq-lbl">Sequence</span><span id="seq-letters"></span></div>

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
  var viewer = $3Dmol.createViewer('viewer', {{backgroundColor:'#0d1117'}});
  viewer.addModel(pdbData, 'pdb');

  var selectedResi = null;
  var labelHandles = [];
  var disulfLines  = [];

  // Precompute feature lookup sets once
  var dSet={{}},oSet={{}},aSet={{}},bSet={{}},gSet={{}},sSet={{}};
  diseaseVars.forEach(function(v)    {{ dSet[v.resi]=1; }});
  otherVars.forEach(function(v)      {{ oSet[v.resi]=1; }});
  activeSites.forEach(function(a)    {{ aSet[a.resi]=1; }});
  bindingSites.forEach(function(b)   {{ bSet[b.resi]=1; }});
  glycanSites.forEach(function(g)    {{ gSet[g.resi]=1; }});
  disulfideBonds.forEach(function(b) {{ sSet[b.pos1]=1; sSet[b.pos2]=1; }});

  // Build sequence strip + cache DOM refs
  var seqEls = {{}};
  (function() {{
    var html='';
    seqData.forEach(function(r) {{
      html += '<span class="seq-aa" id="seq-'+r.resi+'" data-resi="'+r.resi
            + '" title="'+r.resn+' '+r.resi+'">'+r.aa+'</span>';
    }});
    document.getElementById('seq-letters').innerHTML = html;
    seqData.forEach(function(r) {{ seqEls[r.resi]=document.getElementById('seq-'+r.resi); }});
  }})();

  function renderSeq() {{
    var showD=document.getElementById('chk_disease')&&document.getElementById('chk_disease').checked;
    var showO=document.getElementById('chk_other')  &&document.getElementById('chk_other').checked;
    var showA=document.getElementById('chk_active') &&document.getElementById('chk_active').checked;
    var showB=document.getElementById('chk_binding')&&document.getElementById('chk_binding').checked;
    var showG=document.getElementById('chk_glycan') &&document.getElementById('chk_glycan').checked;
    var showS=document.getElementById('chk_disulf') &&document.getElementById('chk_disulf').checked;
    seqData.forEach(function(r) {{
      var el=seqEls[r.resi]; if(!el) return;
      if(r.resi===selectedResi) {{
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
    if(selectedResi!==null&&seqEls[selectedResi])
      seqEls[selectedResi].scrollIntoView({{behavior:'smooth',block:'nearest',inline:'center'}});
  }}

  function renderAll() {{
    viewer.setStyle({{}},{{}});

    if($('#chk_cartoon').prop('checked')) {{
      viewer.setStyle({{}},{{cartoon:{{colorfunc:function(atom){{
        var b=atom.b;
        if(b>=90) return '#0053d6';
        if(b>=70) return '#65cbf3';
        if(b>=50) return '#ffdb13';
        return '#ff7d45';
      }}}}}});
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
    selectedResi=null;
    $('#variant-panel').hide();
    viewer.zoomTo();
    renderAll();
  }};

  function selectVariant(v) {{
    selectedResi=v.resi;
    viewer.zoomTo({{resi:v.resi}},800);
    renderAll();
    showVariantPanel(v);
  }}

  $(document).on('click','.seq-aa',function() {{
    var resi=parseInt($(this).data('resi'));
    selectedResi=resi;
    viewer.zoomTo({{resi:resi}},800);
    var dv=diseaseVars.find(function(v){{ return v.resi===resi; }});
    if(dv) {{ showVariantPanel(dv); }} else {{ $('#variant-panel').hide(); }}
    renderAll();
  }});

  renderAll();
  drawDisulfLines();
  drawLabels();
  viewer.zoomTo();
  viewer.render();

  $('#chk_cartoon').change(function(){{ renderAll(); }});
  $('#chk_disease').change(function(){{ renderAll(); }});
  $('#chk_other').change(function()  {{ renderAll(); }});
  $('#chk_glycan').change(function() {{ renderAll(); }});
  $('#chk_active').change(function() {{ renderAll(); }});
  $('#chk_binding').change(function(){{ renderAll(); }});
  $('#chk_disulf').change(function() {{ drawDisulfLines(); renderAll(); }});
  $('#chk_labels').change(function() {{ drawLabels(); viewer.render(); }});
}});
</script>
</body></html>"""


# ═══════════════════════════════════════════════════════════════════════════════
# Sidebar
# ═══════════════════════════════════════════════════════════════════════════════

with st.sidebar:
    st.markdown("## 🧬 PDB Analyzer")
    st.markdown("---")

    uploaded = st.file_uploader(
        "Upload a PDB file",
        type=["pdb", "ent"],
        help="Standard PDB format files (.pdb or .ent)",
    )

    st.markdown("---")
    st.markdown("### Analysis Settings")
    contact_cutoff = st.slider(
        "Contact cutoff (Å)", min_value=3.0, max_value=7.0,
        value=4.5, step=0.1,
        help="Maximum heavy-atom distance for non-bonded contacts",
    )

    st.markdown("---")
    st.markdown("### BSA Settings")
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

    st.markdown("---")
    st.markdown(
        "**Color scheme**\n"
        "- 🟡 Ligand\n"
        "- 🔵 Contacts\n"
        "- 🟠 H-bonds\n"
        "- 🟣 Pi interactions\n"
        "- 🔵/🔴 Salt bridges"
    )


# ── Setup when PDB uploaded ────────────────────────────────────────────────────
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
            st.image(bf_png, use_container_width=True)
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
        st.caption(
            "BSA(A, B) = ( SASA_A + SASA_B − SASA_AB ) / 2  ·  "
            f"Probe radius: {probe_radius} Å  ·  Sphere points: {n_points_bsa}  ·  "
            "Adjust in sidebar to recompute."
        )
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
                st.image(bsa_png, use_container_width=True)
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
        st.caption(f"Contact cutoff: {contact_cutoff} Å  ·  Adjust in sidebar and results update automatically")
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
                        st.image(fp_bytes, use_container_width=True)
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
                        st.image(d2_bytes, use_container_width=True)
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
            with st.spinner(f"Fetching UniProt annotations for {af_uid}…"):
                uniprot_data = _fetch_uniprot(af_uid)

            all_variants: list[dict] = []
            features: dict = {}
            gene = af_uid
            if uniprot_data:
                all_variants = _parse_variants(uniprot_data)
                features     = _parse_uniprot_features(uniprot_data)
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
