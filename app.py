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

import os
import sys
import tempfile
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
.metric-row { display: flex; gap: 1rem; flex-wrap: wrap; margin-bottom: 1rem; }
.stTabs [data-baseweb="tab-list"] { gap: 6px; }
.stTabs [data-baseweb="tab"] { padding: 8px 18px; border-radius: 6px 6px 0 0; }
</style>
""", unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════════════
# Cached analysis functions (keyed on file bytes so rerun is safe)
# ═══════════════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def _run_structure_analysis(data: bytes) -> dict:
    """Run all summary-level analyses and return plain-Python dicts/lists."""
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        f.write(data)
        p = Path(f.name)
    try:
        ss_structure = ss_load(str(p))
        return {
            "summary":       summarize_structure(p),
            "bfactor":       bfactor_stats(ss_structure),
            "ramachandran":  ramachandran_data(ss_structure),
            "rg":            radius_of_gyration(ss_structure),
            "fasta":         extract_fasta(p),
            "missing":       find_missing_residues(p),
        }
    finally:
        p.unlink(missing_ok=True)


@st.cache_data(show_spinner=False)
def _run_ligand_analysis(data: bytes, cutoff: float):
    """Auto-detect primary ligands and run full interaction analysis."""
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
    """Compute pairwise chain BSA. Returns (results_list, chain_ids_list)."""
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
    """Return all HETATM residues (including excipients) for display."""
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


# ─── Plot helpers (save to temp PNG, return bytes) ────────────────────────────

def _png_bfactor(bf_data: list[dict], title: str) -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_bfactor(bf_data, p, title=title)
    data = p.read_bytes()
    p.unlink(missing_ok=True)
    return data


def _png_ramachandran(rama_data: list[dict], title: str) -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_ramachandran(rama_data, p, title=title)
    data = p.read_bytes()
    p.unlink(missing_ok=True)
    return data


def _png_bsa_matrix(bsa_results: list, chain_ids: list, title: str) -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_bsa_matrix(bsa_results, chain_ids, p, title=title)
    data = p.read_bytes()
    p.unlink(missing_ok=True)
    return data


def _png_fingerprint(report, title: str = "") -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_interaction_summary(report, p)
    data = p.read_bytes()
    p.unlink(missing_ok=True)
    return data


def _png_2d(pdb_path: Path, report, structure) -> bytes:
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        p = Path(f.name)
    plot_ligand_2d(pdb_path, report, p, structure=structure)
    data = p.read_bytes()
    p.unlink(missing_ok=True)
    return data


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

# ── Require upload ─────────────────────────────────────────────────────────────
if not uploaded:
    st.title("🧬 PDB Structure Analyzer")
    st.info(
        "Upload a PDB file in the sidebar to begin.  \n"
        "This app runs full structural analysis including:\n"
        "- Structure metadata, B-factors, Ramachandran, FASTA, Radius of gyration\n"
        "- Ligand interaction analysis (contacts, H-bonds, pi interactions, salt bridges)\n"
        "- Interactive 3D viewer (3Dmol.js) with toggleable interaction layers\n"
        "- 2D LIGPLOT-style diagrams and interaction fingerprints"
    )
    st.stop()

# ── Persist temp file across reruns ───────────────────────────────────────────
pdb_bytes = uploaded.read()

if st.session_state.get("pdb_name") != uploaded.name:
    d = Path(tempfile.mkdtemp())
    p = d / uploaded.name
    p.write_bytes(pdb_bytes)
    st.session_state["pdb_name"] = uploaded.name
    st.session_state["pdb_path"] = str(p)

pdb_path = Path(st.session_state["pdb_path"])

# ── Page title ─────────────────────────────────────────────────────────────────
st.title(f"🧬 {pdb_path.stem.upper()}")
st.caption(f"File: `{uploaded.name}`  ·  {len(pdb_bytes):,} bytes")

# ── Run analyses (cached on file bytes) ────────────────────────────────────────
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

# ═══════════════════════════════════════════════════════════════════════════════
# Top-level metrics bar
# ═══════════════════════════════════════════════════════════════════════════════
c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric("Resolution",
          f"{summary['resolution']} Å" if summary['resolution'] else "N/A")
c2.metric("R-work",
          str(summary["r_work"]) if summary["r_work"] else "N/A")
c3.metric("R-free",
          str(summary["r_free"]) if summary["r_free"] else "N/A")
c4.metric("Residues", f"{summary['n_residues']:,}")
c5.metric("Ligands", summary["n_ligands"])
c6.metric("Rg (Å)", f"{rg:.2f}" if rg else "N/A")

st.divider()

# ═══════════════════════════════════════════════════════════════════════════════
# Tabs
# ═══════════════════════════════════════════════════════════════════════════════
tab_summary, tab_bsa, tab_ligands, tab_viewer, tab_2d = st.tabs([
    "📋 Structure Summary",
    "⬛ Buried Surface Area",
    "🔬 Ligand Analysis",
    "🌐 3D Viewer",
    "🖼 2D Diagrams",
])


# ─── TAB 1: Structure Summary ──────────────────────────────────────────────────
with tab_summary:

    # Metadata table
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

    # Missing residues
    if missing_data:
        st.warning(f"**Missing residue gaps** detected in {len(missing_data)} chain(s)")
        st.dataframe(
            pd.DataFrame(missing_data), use_container_width=True, hide_index=True,
        )
    else:
        st.success("No missing residue gaps detected")

    # All HETATM ligands
    all_ligs = _list_all_ligands(pdb_bytes)
    if all_ligs:
        with st.expander(f"All HETATM residues ({len(all_ligs)})"):
            st.dataframe(
                pd.DataFrame(all_ligs), use_container_width=True, hide_index=True,
            )

    st.divider()

    # ── B-factor ──────────────────────────────────────────────────────────────
    st.subheader("B-factor Profile")
    if bf_data:
        n_high = sum(1 for r in bf_data if r["flag"] == "high")
        n_low  = sum(1 for r in bf_data if r["flag"] == "low")
        st.caption(
            f"{len(bf_data)} Cα residues  ·  "
            f"**{n_high}** high (>2σ)  ·  **{n_low}** low (<2σ)"
        )
        bf_png = _png_bfactor(bf_data, pdb_path.stem.upper())
        st.image(bf_png, use_container_width=True)
        col_dl, _ = st.columns([1, 5])
        with col_dl:
            st.download_button(
                "⬇ B-factor PNG", bf_png,
                file_name=f"{pdb_path.stem}_bfactor.png", mime="image/png",
            )
        with st.expander("B-factor data table"):
            df_bf = pd.DataFrame(bf_data)
            st.dataframe(df_bf, use_container_width=True, hide_index=True)
    else:
        st.info("No B-factor data (no Cα atoms found).")

    st.divider()

    # ── Ramachandran ──────────────────────────────────────────────────────────
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
        rc4.metric("Outliers", n_outlier, delta=f"{pct_ok:.1f}% OK",
                   delta_color="normal")

        rama_png = _png_ramachandran(rama_data, pdb_path.stem.upper())
        col_img, col_dl = st.columns([4, 1])
        with col_img:
            st.image(rama_png)
        with col_dl:
            st.download_button(
                "⬇ Ramachandran PNG", rama_png,
                file_name=f"{pdb_path.stem}_rama.png", mime="image/png",
            )

        if n_outlier:
            with st.expander(f"Outlier residues ({n_outlier})"):
                outliers = [r for r in rama_data if r["region"] == "outlier"]
                st.dataframe(
                    pd.DataFrame(outliers), use_container_width=True, hide_index=True,
                )
    else:
        st.info("No Ramachandran data computed.")

    st.divider()

    # ── FASTA ─────────────────────────────────────────────────────────────────
    st.subheader("FASTA Sequences (from ATOM records)")
    if fasta_str:
        st.code(fasta_str, language=None)
        st.download_button(
            "⬇ Download FASTA", fasta_str,
            file_name=f"{pdb_path.stem}.fasta", mime="text/plain",
        )
    else:
        st.info("No protein sequence extracted.")


# ─── TAB 2: Buried Surface Area ───────────────────────────────────────────────
with tab_bsa:
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
            # Top metrics
            top = bsa_results[0]
            n_iface = sum(1 for r in bsa_results if r["bsa_total"] > 0)
            st.success(
                f"**{n_iface}** interface(s) detected across "
                f"**{len(chain_ids)}** chains  ·  "
                f"Largest: chains {top['chain_1']}–{top['chain_2']} "
                f"({top['bsa_total']:.0f} Å²)"
            )

            # Matrix + bar chart image
            bsa_png = _png_bsa_matrix(
                bsa_results, chain_ids, title=pdb_path.stem.upper()
            )
            st.image(bsa_png, use_container_width=True)
            st.download_button(
                "⬇ BSA matrix PNG", bsa_png,
                file_name=f"{pdb_path.stem}_bsa_matrix.png", mime="image/png",
            )

            st.divider()

            # Full results table
            st.markdown("#### Pairwise BSA table")
            st.caption(
                "`bsa_total` = total interface area (Å²)  ·  "
                "`bsa_on_1/2` = surface buried on each chain individually  ·  "
                "`interface_pct` = % of that chain's SASA buried at this interface"
            )
            import pandas as pd
            df_bsa = pd.DataFrame(bsa_results)
            # Friendly column names
            df_bsa.columns = [
                "Chain 1", "Chain 2",
                "SASA 1 (Å²)", "SASA 2 (Å²)", "SASA complex (Å²)",
                "BSA total (Å²)", "BSA on 1 (Å²)", "BSA on 2 (Å²)",
                "Interface % 1", "Interface % 2",
            ]
            st.dataframe(
                df_bsa.style.background_gradient(
                    subset=["BSA total (Å²)"], cmap="YlOrRd"
                ),
                use_container_width=True, hide_index=True,
            )
            st.download_button(
                "⬇ BSA table CSV",
                df_bsa.to_csv(index=False),
                file_name=f"{pdb_path.stem}_bsa.csv",
                mime="text/csv",
            )

            # Per-chain summary
            st.divider()
            st.markdown("#### Per-chain SASA summary")
            sasa_rows = []
            # Collect unique chain SASAs from results (chain_1 entries)
            seen = {}
            for r in bsa_results:
                if r["chain_1"] not in seen:
                    seen[r["chain_1"]] = r["sasa_1"]
                if r["chain_2"] not in seen:
                    seen[r["chain_2"]] = r["sasa_2"]
            # Also add chains that appear only as chain_2 in no pairs
            for cid in chain_ids:
                if cid not in seen:
                    seen[cid] = None

            for cid in chain_ids:
                sasa_val = seen.get(cid)
                total_buried = sum(
                    r["bsa_on_1"] if r["chain_1"] == cid else r["bsa_on_2"]
                    for r in bsa_results
                    if (r["chain_1"] == cid or r["chain_2"] == cid) and r["bsa_total"] > 0
                )
                sasa_rows.append({
                    "Chain":              cid,
                    "SASA isolated (Å²)": sasa_val,
                    "Total buried (Å²)":  round(total_buried, 1),
                    "% buried":           round(100 * total_buried / sasa_val, 1)
                                          if sasa_val else None,
                })
            st.dataframe(
                pd.DataFrame(sasa_rows), use_container_width=True, hide_index=True,
            )


# ─── TAB 3: Ligand Analysis ────────────────────────────────────────────────────
with tab_ligands:
    st.subheader("Ligand Interaction Analysis")
    st.caption(f"Contact cutoff: {contact_cutoff} Å  ·  "
               "Adjust in sidebar and results update automatically")

    if not reports:
        st.warning(
            "No primary ligands detected.  \n"
            "All HETATM residues are water, crystallographic excipients, "
            "or have fewer than 7 heavy atoms."
        )
    else:
        st.success(f"**{len(reports)}** primary ligand(s) detected and analysed.")

        for report in reports:
            lig_id = (f"{report.ligand_resn} "
                      f"chain {report.ligand_chain} resi {report.ligand_resi}")
            counts = (f"{report.n_contacts} contacts · "
                      f"{report.n_hbonds} H-bonds · "
                      f"{report.n_pi} pi · "
                      f"{report.n_salt_bridges} salt bridges")

            with st.expander(f"**{lig_id}**  —  {counts}", expanded=True):
                it1, it2, it3, it4 = st.tabs([
                    f"Contacts ({report.n_contacts})",
                    f"H-bonds ({report.n_hbonds})",
                    f"Pi interactions ({report.n_pi})",
                    f"Salt bridges ({report.n_salt_bridges})",
                ])

                with it1:
                    if report.contacts:
                        df = pd.DataFrame([c.as_dict() for c in report.contacts])
                        st.dataframe(df, use_container_width=True, hide_index=True)
                        st.download_button(
                            "⬇ CSV", df.to_csv(index=False),
                            file_name=f"{report.ligand_resn}_contacts.csv",
                            mime="text/csv",
                        )
                    else:
                        st.info("No contacts detected.")

                with it2:
                    if report.hbonds:
                        df = pd.DataFrame([h.as_dict() for h in report.hbonds])
                        st.dataframe(df, use_container_width=True, hide_index=True)
                        st.download_button(
                            "⬇ CSV", df.to_csv(index=False),
                            file_name=f"{report.ligand_resn}_hbonds.csv",
                            mime="text/csv",
                        )
                    else:
                        st.info("No hydrogen bonds detected.")

                with it3:
                    if report.pi_interactions:
                        df = pd.DataFrame([p.as_dict() for p in report.pi_interactions])
                        st.dataframe(df, use_container_width=True, hide_index=True)
                        st.download_button(
                            "⬇ CSV", df.to_csv(index=False),
                            file_name=f"{report.ligand_resn}_pi.csv",
                            mime="text/csv",
                        )
                    else:
                        st.info("No pi interactions detected.")

                with it4:
                    if report.salt_bridges:
                        df = pd.DataFrame([s.as_dict() for s in report.salt_bridges])
                        st.dataframe(df, use_container_width=True, hide_index=True)
                        st.download_button(
                            "⬇ CSV", df.to_csv(index=False),
                            file_name=f"{report.ligand_resn}_salt.csv",
                            mime="text/csv",
                        )
                    else:
                        st.info("No salt bridges detected.")


# ─── TAB 3: 3D Viewer ─────────────────────────────────────────────────────────
with tab_viewer:
    st.subheader("Interactive 3D Viewer")

    if not reports:
        st.warning("No primary ligands detected — 3D viewer unavailable.")
    else:
        # Ligand selector
        lig_labels = [
            f"{r.ligand_resn}  chain {r.ligand_chain}:{r.ligand_resi}"
            for r in reports
        ]
        if len(reports) > 1:
            sel_idx = st.selectbox(
                "Select ligand to visualise",
                range(len(reports)),
                format_func=lambda i: lig_labels[i],
            )
        else:
            sel_idx = 0
            st.caption(f"Showing: **{lig_labels[0]}**")

        report = reports[sel_idx]

        st.caption(
            f"{report.summary()}  ·  "
            "Use toggles in the viewer panel to show/hide interaction types. "
            "Drag to rotate · scroll to zoom."
        )

        # Legend row
        lc1, lc2, lc3, lc4, lc5 = st.columns(5)
        lc1.markdown("🟡 **Ligand**")
        lc2.markdown("🔵 **Contacts**")
        lc3.markdown("🟠 **H-bonds**")
        lc4.markdown("🟣 **Pi**")
        lc5.markdown("🔵🔴 **Salt bridges**")

        # Load structure for HTML generation (cached by path string)
        @st.cache_resource
        def _cached_structure(path_str: str):
            return al_load(path_str)

        structure = _cached_structure(str(pdb_path))

        # Render HTML and embed
        html_str = _render_html(
            pdb_path, report, structure,
            width=1100, height=680,
        )
        components.html(html_str, height=720, scrolling=False)

        # Download the HTML file
        st.download_button(
            "⬇ Download standalone HTML viewer",
            html_str,
            file_name=f"{report.ligand_resn}_{report.ligand_chain}{report.ligand_resi}.html",
            mime="text/html",
        )


# ─── TAB 4: 2D Diagrams ───────────────────────────────────────────────────────
with tab_2d:
    st.subheader("2D Interaction Diagrams")

    if not reports:
        st.warning("No primary ligands detected.")
    else:
        @st.cache_resource
        def _cached_structure_2d(path_str: str):
            return al_load(path_str)

        structure_2d = _cached_structure_2d(str(pdb_path))

        for report in reports:
            lig_label = (
                f"{report.ligand_resn} "
                f"chain {report.ligand_chain}:{report.ligand_resi}"
            )
            st.markdown(f"### {lig_label}")
            st.caption(report.summary())

            col_fp, col_2d = st.columns(2, gap="large")

            with col_fp:
                st.markdown("#### Interaction Fingerprint")
                st.caption(
                    "Dot size = contact count · Dot shade = proximity · "
                    "Number = count when >1"
                )
                if (report.n_contacts + report.n_hbonds +
                        report.n_pi + report.n_salt_bridges) > 0:
                    fp_bytes = _png_fingerprint(report)
                    st.image(fp_bytes, use_container_width=True)
                    st.download_button(
                        "⬇ Fingerprint PNG", fp_bytes,
                        file_name=f"{report.ligand_resn}_fingerprint.png",
                        mime="image/png",
                        key=f"dl_fp_{report.ligand_resn}_{report.ligand_resi}",
                    )
                else:
                    st.info("No interactions to plot.")

            with col_2d:
                st.markdown("#### LIGPLOT-style 2D Diagram")
                st.caption(
                    "Ligand in centre (CPK colours) · "
                    "Protein residues at periphery · "
                    "Spoked arcs = hydrophobic contacts"
                )
                if (report.n_contacts + report.n_hbonds +
                        report.n_pi + report.n_salt_bridges) > 0:
                    d2_bytes = _png_2d(pdb_path, report, structure_2d)
                    st.image(d2_bytes, use_container_width=True)
                    st.download_button(
                        "⬇ 2D diagram PNG", d2_bytes,
                        file_name=f"{report.ligand_resn}_2d.png",
                        mime="image/png",
                        key=f"dl_2d_{report.ligand_resn}_{report.ligand_resi}",
                    )
                else:
                    st.info("No interactions to plot.")

            st.divider()
