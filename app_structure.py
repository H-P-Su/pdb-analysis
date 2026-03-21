#!/usr/bin/env python3
"""
app_structure.py — PDB Structure Analyzer.

Integrates:
  • summarize_structures.py  — metadata, B-factors, Ramachandran, FASTA, Rg
  • analyze_ligands.py       — contacts, H-bonds, pi interactions, salt bridges
  • visualize_interactions.py — 3D HTML viewer (3Dmol.js) + 2D LIGPLOT diagrams

Run:
    streamlit run app_structure.py
"""

from __future__ import annotations

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

st.markdown("""
<style>
.stTabs [data-baseweb="tab-list"] { gap: 6px; }
.stTabs [data-baseweb="tab"] { padding: 8px 18px; border-radius: 6px 6px 0 0; }
</style>
""", unsafe_allow_html=True)

_NO_PDB = "Upload a PDB file above to get started."


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


# ── Plot helpers ───────────────────────────────────────────────────────────────

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

    st.title(f"🔬 {pdb_path.stem.upper()}")
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
# Tabs
# ═══════════════════════════════════════════════════════════════════════════════
tab_summary, tab_bsa, tab_viewer, tab_2d = st.tabs([
    "📋 Structure Summary",
    "⬛ Buried Surface Area",
    "🌐 3D Viewer",
    "🖼 2D Diagrams",
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
                        "Chain": cid,
                        "SASA isolated (A\u00b2)": round(sv, 1) if sv else None,
                        "Total buried (A\u00b2)": round(tb, 1),
                        "% buried": round(100 * tb / sv, 1) if sv else None,
                    })
                st.dataframe(pd.DataFrame(sasa_rows), hide_index=True)


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
