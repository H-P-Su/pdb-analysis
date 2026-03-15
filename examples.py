#!/usr/bin/env python3
"""
examples.py — Demonstration of all project functionality.

Downloads four demo structures and runs every analysis feature,
saving all outputs to demo_output/.

Run with:
    python3 examples.py

Structures used:
    4HHB  — Oxy-hemoglobin  (chains A–D, HEM ligands)
    2HHB  — Deoxy-hemoglobin (good RMSD pair with 4HHB)
    6D1Y  — Kinase with small-molecule inhibitor FQJ
    6D1Z  — Same kinase, different conformation (RMSD reference)
    6D20  — Same kinase series
    1D3Z  — Ubiquitin NMR ensemble (10 models)
"""

import sys
import traceback
from pathlib import Path

# ─── Helpers ──────────────────────────────────────────────────────────────────

DEMO_DIR  = Path("demo_output")
FILES_DIR = Path("Files")
DEMO_DIR.mkdir(exist_ok=True)
FILES_DIR.mkdir(exist_ok=True)

_section_num = [0]

def banner(title: str):
    _section_num[0] += 1
    n = _section_num[0]
    bar = "═" * 62
    print(f"\n{bar}")
    print(f"  {n}. {title}")
    print(f"{bar}\n")

def sub(title: str):
    print(f"\n  ── {title}")

def run(fn, *args, label="", **kwargs):
    """Call fn(*args, **kwargs), print label, catch and report any errors."""
    tag = label or fn.__name__
    try:
        result = fn(*args, **kwargs)
        return result
    except Exception:
        print(f"\n  [ERROR in {tag}]")
        traceback.print_exc()
        return None


# ─── 1. DOWNLOADS ─────────────────────────────────────────────────────────────

banner("DOWNLOADS")

from download_pdb import download_entry

DEMO_IDS = {
    "4HHB": "Oxy-hemoglobin with HEM ligands (chains A–D)",
    "2HHB": "Deoxy-hemoglobin — RMSD pair for 4HHB",
    "6D1Y": "Kinase with small-molecule inhibitor FQJ",
    "6D1Z": "Same kinase, apo/alt conformation",
    "6D20": "Same kinase series",
    "1D3Z": "Ubiquitin NMR ensemble (10 models)",
}

for pdb_id, desc in DEMO_IDS.items():
    print(f"  {pdb_id}: {desc}")
    download_entry(pdb_id, FILES_DIR,
                   fmt_pdb=True, fmt_cif=False, fmt_seq=True, overwrite=False)

print("\n  All demo structures ready in Files/")


# ─── 2. STRUCTURE SUMMARIES ───────────────────────────────────────────────────

banner("STRUCTURE SUMMARIES")

from summarize_structures import (
    summarize_structure, batch_summary, find_missing_residues,
    extract_fasta, bfactor_stats, radius_of_gyration,
    ramachandran_data, plot_bfactor, plot_ramachandran,
    load_structure, print_summary,
)

# 2a. Single structure summary
sub("Single structure summary: 4HHB (hemoglobin)")
row = summarize_structure(FILES_DIR / "4HHB.pdb")
print_summary(row)

# 2b. Batch summary → CSV
sub("Batch summary → demo_output/batch_summary.csv")
all_paths = [FILES_DIR / f"{pid}.pdb" for pid in DEMO_IDS]
batch_summary(all_paths, DEMO_DIR / "batch_summary.csv")

# 2c. Missing residues
sub("Missing residues: 6D1Y")
missing = find_missing_residues(FILES_DIR / "6D1Y.pdb")
if missing:
    for m in missing:
        print(f"    Chain {m['chain']}: {m['observed_residues']}/{m['expected_residues']} "
              f"modeled ({m['pct_modeled']}%,  {m['missing_count']} missing)")
else:
    print("    No missing residues detected.")

# 2d. FASTA from coordinates
sub("FASTA from ATOM records: 4HHB chain A (first 3 lines)")
fasta = extract_fasta(FILES_DIR / "4HHB.pdb", chain_id="A")
for line in fasta.splitlines()[:3]:
    print(f"    {line}")

# 2e. Radius of gyration
sub("Radius of gyration")
for pid in ["4HHB", "6D1Y", "1D3Z"]:
    s = load_structure(FILES_DIR / f"{pid}.pdb")
    rg = radius_of_gyration(s)
    print(f"    {pid}  Rg = {rg:.2f} Å")


# ─── 3. GEOMETRY ANALYSIS ─────────────────────────────────────────────────────

banner("GEOMETRY ANALYSIS")

# 3a. B-factor
sub("B-factor profile: 6D1Y → demo_output/6D1Y_bfactor.png")
s_6d1y = load_structure(FILES_DIR / "6D1Y.pdb")
bf_data = bfactor_stats(s_6d1y)
high = [r for r in bf_data if r["flag"] == "high"]
low  = [r for r in bf_data if r["flag"] == "low"]
print(f"    {len(bf_data)} residues  |  {len(high)} high-B  |  {len(low)} low-B")
if high:
    print(f"    High-B residues: " +
          ", ".join(f"{r['resn']}{r['resi']}({r['bfactor']:.0f})" for r in high[:5]))
plot_bfactor(bf_data, DEMO_DIR / "6D1Y_bfactor.png", title="6D1Y")

# 3b. Ramachandran
sub("Ramachandran plot: 6D1Y → demo_output/6D1Y_rama.png")
rama = ramachandran_data(s_6d1y)
by_region = {}
for r in rama:
    by_region.setdefault(r["region"], []).append(r)
for region, rows in sorted(by_region.items()):
    pct = 100 * len(rows) / len(rama)
    print(f"    {region:10s}: {len(rows):4d} residues  ({pct:.1f}%)")
plot_ramachandran(rama, DEMO_DIR / "6D1Y_rama.png", title="6D1Y")

# 3c. Ramachandran for 4HHB to compare
sub("Ramachandran plot: 4HHB → demo_output/4HHB_rama.png")
s_4hhb = load_structure(FILES_DIR / "4HHB.pdb")
rama_4hhb = ramachandran_data(s_4hhb)
outliers_4hhb = [r for r in rama_4hhb if r["region"] == "outlier"]
pct_ok = 100 * (1 - len(outliers_4hhb) / len(rama_4hhb)) if rama_4hhb else 0
print(f"    {len(rama_4hhb)} residues  |  {len(outliers_4hhb)} outliers  |  {pct_ok:.1f}% allowed")
plot_ramachandran(rama_4hhb, DEMO_DIR / "4HHB_rama.png", title="4HHB")


# ─── 4. LIGAND ANALYSIS ───────────────────────────────────────────────────────

banner("LIGAND ANALYSIS  (analyze_ligands.py)")

from analyze_ligands import (
    load_structure as al_load, get_ligands,
    find_contacts, find_hydrogen_bonds,
    find_pi_interactions, find_salt_bridges,
    analyze_ligand, to_csv, print_table,
)

s_6d1y_al = al_load(FILES_DIR / "6D1Y.pdb")

# 4a. List ligands
sub("List ligands: 6D1Y")
ligs = get_ligands(s_6d1y_al)
for res in ligs:
    atoms = sum(1 for _ in res.get_atoms())
    print(f"    {res.get_resname():<6} chain {res.get_parent().id}  "
          f"resi {res.id[1]}  ({atoms} atoms)")

# 4b. Contacts
sub("Contacts: FQJ ↔ polymer (cutoff 4.5 Å)")
contacts = find_contacts(s_6d1y_al, "resn FQJ", "polymer", cutoff=4.5)
print(f"    {len(contacts)} contacts found  (closest: "
      f"{contacts[0].resn2}{contacts[0].resi2} {contacts[0].distance} Å)")
to_csv(contacts, DEMO_DIR / "6D1Y_FQJ_contacts.csv")

# 4c. H-bonds
sub("Hydrogen bonds: polymer ↔ FQJ")
hbonds_fwd = find_hydrogen_bonds(s_6d1y_al, "resn FQJ", "polymer")
hbonds_rev = find_hydrogen_bonds(s_6d1y_al, "polymer",  "resn FQJ")
hbonds = hbonds_fwd + hbonds_rev
print(f"    {len(hbonds)} H-bonds found")
print_table(hbonds[:5])
to_csv(hbonds, DEMO_DIR / "6D1Y_FQJ_hbonds.csv")

# 4d. Pi interactions
sub("Pi interactions: polymer ↔ FQJ")
pi = find_pi_interactions(s_6d1y_al, "polymer", "resn FQJ",
                           pdb_path=FILES_DIR / "6D1Y.pdb")
print(f"    {len(pi)} pi interactions found")
print_table(pi)
to_csv(pi, DEMO_DIR / "6D1Y_FQJ_pi.csv")

# 4e. Salt bridges (protein–protein)
sub("Salt bridges: polymer ↔ polymer (4HHB)")
s_4hhb_al = al_load(FILES_DIR / "4HHB.pdb")
sb = find_salt_bridges(s_4hhb_al, "polymer", "polymer")
print(f"    {len(sb)} salt bridges found  (shortest: "
      f"{sb[0].cation_resn}{sb[0].cation_resi} ↔ "
      f"{sb[0].anion_resn}{sb[0].anion_resi}  {sb[0].distance} Å)")
to_csv(sb, DEMO_DIR / "4HHB_salt_bridges.csv")

# 4f. Full ligand report: HEM in 4HHB
sub("Full ligand report: HEM in 4HHB (chain A)")
reports_hhb = analyze_ligand(s_4hhb_al, "resn HEM and chain A",
                              protein_sel="polymer and chain A")
for report in reports_hhb:
    print(f"    {report.summary()}")

# 4g. Full ligand report: FQJ in 6D1Y
sub("Full ligand report: FQJ in 6D1Y")
reports_6d1y = analyze_ligand(s_6d1y_al, "resn FQJ",
                               pdb_path=FILES_DIR / "6D1Y.pdb")
for report in reports_6d1y:
    print(f"    {report.summary()}")


# ─── 5. VISUALIZATION ─────────────────────────────────────────────────────────

banner("VISUALIZATION  (visualize_interactions.py)")

from visualize_interactions import visualize

# 5a. Kinase inhibitor: interactive HTML + fingerprint PNG
sub("Kinase inhibitor FQJ: 3D HTML + fingerprint PNG")
results = visualize(
    FILES_DIR / "6D1Y.pdb",
    ligand_sel  = "resn FQJ",
    protein_sel = "polymer",
    outdir      = DEMO_DIR,
    prefix      = "demo_",
    width=900, height=700,
)
for r in results:
    print(f"    HTML: {r['html']}")
    print(f"    PNG:  {r['png']}")

# 5b. Hemoglobin HEM ligand (chain A only)
sub("Hemoglobin HEM (chain A): 3D HTML + fingerprint PNG")
results_hhb = visualize(
    FILES_DIR / "4HHB.pdb",
    ligand_sel  = "resn HEM and chain A",
    protein_sel = "polymer and chain A",
    outdir      = DEMO_DIR,
    prefix      = "demo_",
    width=900, height=700,
)
for r in results_hhb:
    print(f"    HTML: {r['html']}")
    print(f"    PNG:  {r['png']}")


# ─── 6. RMSD ANALYSIS ────────────────────────────────────────────────────────

banner("RMSD ANALYSIS  (analyze_rmsd.py)")

from analyze_rmsd import (
    calculate_rmsd, per_residue_rmsd, superpose,
    pairwise_rmsd_matrix, nmr_ensemble_rmsd,
    plot_rmsd_matrix, plot_per_residue_rmsd, plot_ensemble_rmsd,
    save_csv,
)

# 6a. Simple RMSD between two structures
sub("RMSD: 4HHB vs 2HHB (oxy vs deoxy hemoglobin, chain A)")
rmsd = calculate_rmsd(FILES_DIR / "4HHB.pdb", FILES_DIR / "2HHB.pdb",
                       chain1="A", chain2="A")
print(f"    Cα RMSD (chain A): {rmsd:.3f} Å")

rmsd_full = calculate_rmsd(FILES_DIR / "4HHB.pdb", FILES_DIR / "2HHB.pdb")
print(f"    Cα RMSD (all chains): {rmsd_full:.3f} Å")

# 6b. RMSD by residue number (same numbering scheme)
sub("RMSD: 6D1Z vs 6D20 (two kinase conformers)")
rmsd_close = calculate_rmsd(FILES_DIR / "6D1Z.pdb", FILES_DIR / "6D20.pdb")
print(f"    Cα RMSD: {rmsd_close:.3f} Å  (very similar conformations expected)")

# 6c. Per-residue RMSD + plot
sub("Per-residue RMSD: 4HHB vs 2HHB chain A → demo_output/hhb_per_residue.png")
pr = per_residue_rmsd(FILES_DIR / "4HHB.pdb", FILES_DIR / "2HHB.pdb",
                       chain1="A", chain2="A")
if pr:
    max_r = max(pr, key=lambda r: r["rmsd_ca"])
    mean_r = sum(r["rmsd_ca"] for r in pr) / len(pr)
    print(f"    {len(pr)} aligned residues  |  mean {mean_r:.2f} Å  |  "
          f"max {max_r['rmsd_ca']:.2f} Å at {max_r['resn']}{max_r['resi']}")
    save_csv(pr, DEMO_DIR / "hhb_per_residue_rmsd.csv")
    plot_per_residue_rmsd(pr, DEMO_DIR / "hhb_per_residue_rmsd.png",
                          title="4HHB vs 2HHB (chain A)", threshold=1.5)

# 6d. Superposition: save transformed structure
sub("Superpose 2HHB onto 4HHB (chain A) → demo_output/2HHB_superposed.pdb")
rmsd_sup, _ = superpose(
    FILES_DIR / "4HHB.pdb", FILES_DIR / "2HHB.pdb",
    chain1="A", chain2="A",
    output_path=DEMO_DIR / "2HHB_superposed.pdb",
)
print(f"    RMSD after superposition: {rmsd_sup:.3f} Å")

# 6e. Pairwise RMSD matrix → heatmap
sub("Pairwise RMSD matrix: kinase series → demo_output/kinase_rmsd_matrix.png")
kinase_paths = [FILES_DIR / f"{pid}.pdb" for pid in ["6D1Y", "6D1Z", "6D20"]]
matrix, labels = pairwise_rmsd_matrix(kinase_paths)
plot_rmsd_matrix(matrix, labels, DEMO_DIR / "kinase_rmsd_matrix.png",
                 title="Kinase Series — Pairwise Cα RMSD (Å)")

# 6f. NMR ensemble RMSD
sub("NMR ensemble RMSD: 1D3Z (ubiquitin) → demo_output/1D3Z_ensemble.png")
ensemble = nmr_ensemble_rmsd(FILES_DIR / "1D3Z.pdb")
if ensemble:
    non_ref = [r for r in ensemble if r["rmsd"] > 0]
    mean_e = sum(r["rmsd"] for r in non_ref) / len(non_ref) if non_ref else 0
    print(f"    {len(ensemble)} models  |  mean RMSD vs model 0: {mean_e:.3f} Å")
    save_csv(ensemble, DEMO_DIR / "1D3Z_ensemble_rmsd.csv")
    plot_ensemble_rmsd(ensemble, DEMO_DIR / "1D3Z_ensemble_rmsd.png",
                       title="1D3Z Ubiquitin NMR Ensemble")


# ─── Summary ──────────────────────────────────────────────────────────────────

banner("DEMO COMPLETE")

outputs = sorted(DEMO_DIR.iterdir())
print(f"  {len(outputs)} files written to {DEMO_DIR}/\n")
for f in outputs:
    size = f.stat().st_size
    size_str = f"{size/1024:.0f} KB" if size > 1024 else f"{size} B"
    print(f"    {f.name:<45}  {size_str:>8}")

print("""
  ┌─────────────────────────────────────────────────────────┐
  │  Scripts available:                                     │
  │    download_pdb.py          Download PDB / CIF / FASTA  │
  │    analyze_ligands.py       Contacts, H-bonds, Pi, Salt │
  │    visualize_interactions.py  HTML 3D viewer + PNG      │
  │    summarize_structures.py  Metadata, B-factor, Rama    │
  │    analyze_rmsd.py          RMSD, superposition, NMR    │
  └─────────────────────────────────────────────────────────┘
""")
