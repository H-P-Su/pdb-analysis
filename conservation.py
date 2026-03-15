#!/usr/bin/env python3
"""
conservation.py — Sequence conservation mapping from a user-supplied FASTA file.

Aligns each sequence in the FASTA to the reference protein sequence extracted from
the PDB, then computes per-residue conservation scores and maps them onto the
structure (B-factor column).  No BLAST or network access required.

Algorithm
---------
1. Extract the reference sequence from ATOM records (one letter per residue).
2. Read all sequences from the input FASTA file.
3. Pairwise-align each query to the reference using BioPython PairwiseAligner
   (BLOSUM62, global alignment, gap penalties matched to BLAST defaults).
4. Build a position-specific alignment matrix (one column per reference position).
5. Compute per-column Shannon entropy → conservation score = 1 − H/H_max.
   Score 1.0 = fully conserved, 0.0 = fully variable.
6. Output:
     <stem>_conservation.png        bar chart of per-residue scores
     <stem>_conservation.pdb        PDB with conservation in B-factor column
     <stem>_conservation.csv        residue, aa, score, entropy, gap_fraction

Functions
---------
read_fasta(path)                    → list of (header, seq) tuples
extract_ref_sequence(pdb_path, chain)
                                    → list of (resnum, aa) tuples from ATOM records
align_sequences(ref_seq, queries)   → MSA-like column list per ref position
conservation_scores(columns)        → np.ndarray of scores (0–1) per position
write_conservation_pdb(pdb_path, resnum_scores, out_path)
plot_conservation(resnum_scores, out_path, title)
conservation_report(pdb_path, fasta_path, chain, out_dir)

CLI
---
python3 conservation.py 1abc.pdb sequences.fasta
python3 conservation.py 1abc.pdb sequences.fasta --chain A
python3 conservation.py 1abc.pdb sequences.fasta --out results/
python3 conservation.py 1abc.pdb sequences.fasta --plot --no-pdb
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from Bio import Align
from Bio.PDB import PDBParser, MMCIFParser

# ─── Constants ────────────────────────────────────────────────────────────────

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "SEC": "U",
}

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
GAP = "-"

# ─── FASTA parsing ────────────────────────────────────────────────────────────

def read_fasta(path: Path) -> list[tuple[str, str]]:
    """
    Parse a multi-FASTA file.

    Returns a list of (header, sequence) tuples. Sequences are uppercased and
    whitespace-stripped; lines starting with ';' are treated as comments.
    """
    records: list[tuple[str, str]] = []
    header = ""
    seq_parts: list[str] = []

    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq_parts).upper()))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line.replace(" ", ""))

    if header:
        records.append((header, "".join(seq_parts).upper()))

    return records


# ─── Reference sequence from PDB ─────────────────────────────────────────────

def extract_ref_sequence(
    pdb_path: Path,
    chain: Optional[str] = None,
) -> list[tuple[int, str]]:
    """
    Extract the protein sequence from ATOM records of a PDB/mmCIF file.

    Returns a list of (residue_number, one_letter_code) tuples, one per unique
    residue position.  If chain is None, uses the first chain found.
    Non-standard residues are mapped via THREE_TO_ONE; unknowns become 'X'.
    """
    path = Path(pdb_path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("ref", str(path))
    model     = next(iter(structure))

    if chain is None:
        # Pick first chain with protein residues
        for ch in model:
            has_protein = any(
                res.get_resname() in THREE_TO_ONE
                for res in ch
                if res.id[0] == " "
            )
            if has_protein:
                chain = ch.id
                break
        if chain is None:
            raise ValueError(f"No protein chain found in {path.name}")

    if chain not in model:
        raise ValueError(f"Chain '{chain}' not found in {path.name}. "
                         f"Available: {[c.id for c in model]}")

    seen: set[int] = set()
    residues: list[tuple[int, str]] = []
    for res in model[chain]:
        if res.id[0] != " ":   # skip HETATM and water
            continue
        resnum = res.id[1]
        if resnum in seen:
            continue
        seen.add(resnum)
        aa = THREE_TO_ONE.get(res.get_resname(), "X")
        residues.append((resnum, aa))

    if not residues:
        raise ValueError(f"No ATOM residues found in chain {chain} of {path.name}")

    return residues


# ─── Pairwise alignment → column matrix ───────────────────────────────────────

def _make_aligner() -> Align.PairwiseAligner:
    """BioPython PairwiseAligner configured for protein global alignment."""
    aligner = Align.PairwiseAligner()
    aligner.mode              = "global"
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score    = -11.0   # match BLAST protein defaults
    aligner.extend_gap_score  = -1.0
    return aligner


def align_sequences(
    ref_residues: list[tuple[int, str]],
    queries: list[tuple[str, str]],
) -> list[list[str]]:
    """
    Align each query sequence to the reference using BLOSUM62 global alignment.

    Returns a list of columns, one per reference residue position.
    Each column is a list of amino-acid characters (or '-' for gaps) from all
    query sequences that aligned to that reference position.

    Parameters
    ----------
    ref_residues : list of (resnum, aa) from extract_ref_sequence()
    queries      : list of (header, sequence) from read_fasta()
    """
    ref_seq = "".join(aa for _, aa in ref_residues)
    n_ref   = len(ref_seq)

    # columns[i] = list of query residues aligned to ref position i
    columns: list[list[str]] = [[] for _ in range(n_ref)]

    aligner = _make_aligner()

    for header, qseq in queries:
        # Keep only standard amino acid characters
        qseq_clean = "".join(c if c in AMINO_ACIDS else "X" for c in qseq)
        if not qseq_clean:
            print(f"  [warn] Skipping empty sequence: {header[:60]}", file=sys.stderr)
            continue

        try:
            alignments = aligner.align(ref_seq, qseq_clean)
            aln        = next(iter(alignments))  # best alignment
        except StopIteration:
            print(f"  [warn] No alignment for: {header[:60]}", file=sys.stderr)
            continue

        # aln.aligned: tuple of (ref_blocks, query_blocks)
        # Build per-ref-position character by walking the alignment
        ref_aln, qry_aln = _alignment_to_strings(aln, ref_seq, qseq_clean)

        ref_pos = 0   # position in ref_seq (0-indexed)
        for r_char, q_char in zip(ref_aln, qry_aln):
            if r_char != GAP:
                if ref_pos < n_ref:
                    columns[ref_pos].append(q_char)
                ref_pos += 1

    return columns


def _alignment_to_strings(aln, ref_seq: str, qry_seq: str) -> tuple[str, str]:
    """
    Convert a BioPython alignment object to two gapped strings.
    Works with both older (format()) and newer (aligned) APIs.
    """
    try:
        # BioPython ≥ 1.80: aln.aligned gives block coordinates
        ref_blocks, qry_blocks = aln.aligned
        ref_str, qry_str = [], []
        ref_pos, qry_pos = 0, 0

        for (r_start, r_end), (q_start, q_end) in zip(ref_blocks, qry_blocks):
            # gap in ref before this block
            if r_start > ref_pos:
                ref_str.append(GAP * (r_start - ref_pos))
                qry_str.append(qry_seq[qry_pos:q_start])
            ref_str.append(ref_seq[r_start:r_end])
            qry_str.append(qry_seq[q_start:q_end])
            ref_pos, qry_pos = r_end, q_end

        # trailing gaps
        if ref_pos < len(ref_seq):
            ref_str.append(ref_seq[ref_pos:])
            qry_str.append(GAP * (len(ref_seq) - ref_pos))

        return "".join(ref_str), "".join(qry_str)

    except AttributeError:
        # Fallback: use string representation
        lines  = str(aln).splitlines()
        r_line = lines[0] if lines else ""
        q_line = lines[2] if len(lines) > 2 else ""
        return r_line, q_line


# ─── Conservation scoring ─────────────────────────────────────────────────────

def conservation_scores(
    columns: list[list[str]],
    include_gaps_in_entropy: bool = False,
) -> np.ndarray:
    """
    Compute per-position conservation score from alignment columns.

    Score = 1 - H / H_max   where H is Shannon entropy and H_max = log2(20).
    Score 1.0 = fully conserved, 0.0 = maximum variability.
    Positions with no aligned sequences (or all gaps) score 0.

    Parameters
    ----------
    columns : list of columns from align_sequences()
    include_gaps_in_entropy : if True, gap '-' counts as a 21st symbol
    """
    H_MAX = math.log2(20)
    scores = np.zeros(len(columns))

    for i, col in enumerate(columns):
        if not col:
            continue
        symbols = col if include_gaps_in_entropy else [c for c in col if c != GAP]
        if not symbols:
            continue
        counts = Counter(symbols)
        total  = sum(counts.values())
        H      = -sum((c / total) * math.log2(c / total) for c in counts.values())
        scores[i] = max(0.0, 1.0 - H / H_MAX)

    return scores


def gap_fractions(columns: list[list[str]]) -> np.ndarray:
    """Per-position fraction of aligned sequences that have a gap."""
    fracs = np.zeros(len(columns))
    for i, col in enumerate(columns):
        if col:
            fracs[i] = sum(1 for c in col if c == GAP) / len(col)
    return fracs


# ─── Output: CSV ──────────────────────────────────────────────────────────────

def write_csv(
    ref_residues: list[tuple[int, str]],
    scores: np.ndarray,
    columns: list[list[str]],
    out_path: Path,
) -> None:
    """Write per-residue conservation table as CSV."""
    H_MAX = math.log2(20)
    fracs = gap_fractions(columns)

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["resnum", "aa", "conservation", "entropy", "gap_fraction",
                         "n_sequences"])
        for i, (resnum, aa) in enumerate(ref_residues):
            col  = columns[i]
            syms = [c for c in col if c != GAP]
            if syms:
                counts = Counter(syms)
                total  = sum(counts.values())
                H      = -sum((c / total) * math.log2(c / total)
                              for c in counts.values())
            else:
                H = 0.0
            writer.writerow([
                resnum, aa,
                f"{scores[i]:.4f}",
                f"{H:.4f}",
                f"{fracs[i]:.4f}",
                len(col),
            ])


# ─── Output: B-factor PDB ─────────────────────────────────────────────────────

def write_conservation_pdb(
    pdb_path: Path,
    resnum_scores: dict[int, float],
    out_path: Path,
    chain: Optional[str] = None,
) -> None:
    """
    Write a copy of the PDB with conservation scores in the B-factor column.

    Score is scaled to 0–100 (so PyMOL's default B-factor colour ramp applies).
    Residues not in resnum_scores keep their original B-factor.
    If chain is specified, only that chain is modified.
    """
    with open(pdb_path) as fi, open(out_path, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try:
                    line_chain = line[21]
                    resi       = int(line[22:26])
                    if (chain is None or line_chain == chain) and resi in resnum_scores:
                        val  = resnum_scores[resi] * 100.0   # 0–100
                        line = line[:60] + f"{val:6.2f}" + line[66:]
                except (ValueError, IndexError):
                    pass
            fo.write(line)


# ─── Output: plot ─────────────────────────────────────────────────────────────

def plot_conservation(
    ref_residues: list[tuple[int, str]],
    scores: np.ndarray,
    out_path: Path,
    title: str = "Sequence conservation",
    n_seqs: int = 0,
) -> None:
    """
    Bar chart of per-residue conservation scores, coloured from red (variable)
    to blue (conserved).  Saves as PNG.
    """
    resnums = [r for r, _ in ref_residues]
    cmap    = plt.cm.RdYlBu  # red=0 (variable) → blue=1 (conserved)
    colors  = [cmap(s) for s in scores]

    fig, ax = plt.subplots(figsize=(max(10, len(resnums) / 8), 4))
    bars = ax.bar(resnums, scores, color=colors, width=1.0, linewidth=0)

    ax.set_xlabel("Residue number")
    ax.set_ylabel("Conservation score")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"{title}" + (f"  (n={n_seqs} sequences)" if n_seqs else ""))

    # Colorbar legend
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.01)
    cbar.set_label("Score")
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(["Variable", "0.5", "Conserved"])

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ─── High-level entry point ───────────────────────────────────────────────────

def conservation_report(
    pdb_path: Path,
    fasta_path: Path,
    chain: Optional[str] = None,
    out_dir: Optional[Path] = None,
    write_pdb: bool = True,
    write_plot: bool = True,
) -> dict:
    """
    Full conservation analysis pipeline.

    Returns a dict with keys:
      ref_residues  — list of (resnum, aa)
      scores        — np.ndarray of conservation scores
      columns       — list of alignment columns
      n_seqs        — number of query sequences used
      out_csv       — Path to CSV output
      out_pdb       — Path to B-factor PDB (or None)
      out_plot      — Path to PNG plot (or None)
    """
    pdb_path  = Path(pdb_path)
    fasta_path = Path(fasta_path)
    stem      = pdb_path.stem

    if out_dir is None:
        out_dir = pdb_path.parent
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Reference sequence
    print(f"  Extracting reference sequence from {pdb_path.name} …")
    ref_residues = extract_ref_sequence(pdb_path, chain=chain)
    ref_seq      = "".join(aa for _, aa in ref_residues)
    used_chain   = chain or "auto"
    print(f"  Reference: {len(ref_residues)} residues  (chain {used_chain})")
    print(f"  Sequence:  {ref_seq[:60]}{'…' if len(ref_seq) > 60 else ''}")

    # 2. Query sequences
    print(f"  Reading FASTA: {fasta_path.name} …")
    queries = read_fasta(fasta_path)
    if not queries:
        raise ValueError(f"No sequences found in {fasta_path}")
    print(f"  Found {len(queries)} query sequences")

    # 3. Align
    print("  Aligning sequences …")
    columns = align_sequences(ref_residues, queries)
    n_aligned = sum(1 for col in columns if col)
    print(f"  Aligned {n_aligned}/{len(ref_residues)} reference positions covered")

    # 4. Score
    scores = conservation_scores(columns)
    fracs  = gap_fractions(columns)
    mean_s = scores.mean()
    n_high = int((scores >= 0.8).sum())
    n_low  = int((scores <= 0.2).sum())
    print(f"  Conservation: mean={mean_s:.3f}  "
          f"highly conserved (≥0.8): {n_high}  "
          f"variable (≤0.2): {n_low}")

    resnum_scores = {resnum: scores[i] for i, (resnum, _) in enumerate(ref_residues)}

    # 5. CSV
    out_csv = out_dir / f"{stem}_conservation.csv"
    write_csv(ref_residues, scores, columns, out_csv)
    print(f"  CSV → {out_csv.name}")

    # 6. B-factor PDB
    out_pdb = None
    if write_pdb:
        out_pdb = out_dir / f"{stem}_conservation.pdb"
        write_conservation_pdb(pdb_path, resnum_scores, out_pdb, chain=chain)
        print(f"  PDB → {out_pdb.name}")
        print(f"  PyMOL:  load {out_pdb.name};  spectrum b, red_white_blue, minimum=0, maximum=100")

    # 7. Plot
    out_plot = None
    if write_plot:
        out_plot = out_dir / f"{stem}_conservation.png"
        title    = f"Sequence conservation — {stem}"
        plot_conservation(ref_residues, scores, out_plot,
                          title=title, n_seqs=len(queries))
        print(f"  Plot → {out_plot.name}")

    return {
        "ref_residues": ref_residues,
        "scores":       scores,
        "columns":      columns,
        "n_seqs":       len(queries),
        "out_csv":      out_csv,
        "out_pdb":      out_pdb,
        "out_plot":     out_plot,
    }


# ─── CLI ──────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="conservation.py",
        description=(
            "Map sequence conservation onto a PDB structure using a user-supplied "
            "FASTA file.  No BLAST or network access required.\n\n"
            "Each sequence in the FASTA is aligned to the reference protein sequence "
            "extracted from the PDB using BLOSUM62 global alignment (BioPython). "
            "Per-residue Shannon entropy is computed and converted to a conservation "
            "score (1 = fully conserved, 0 = fully variable)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 conservation.py 1abc.pdb homologs.fasta\n"
            "  python3 conservation.py 1abc.pdb homologs.fasta --chain A\n"
            "  python3 conservation.py 1abc.pdb homologs.fasta --out results/\n"
            "  python3 conservation.py 1abc.pdb homologs.fasta --no-pdb\n"
        ),
    )
    p.add_argument("pdb",   type=Path, help="Input PDB or mmCIF file")
    p.add_argument("fasta", type=Path, help="Multi-FASTA file of homologous sequences")
    p.add_argument("--chain",  default=None,
                   help="Chain to use as reference (default: first protein chain)")
    p.add_argument("--out",    type=Path, default=None,
                   help="Output directory (default: same as PDB file)")
    p.add_argument("--no-pdb",  action="store_true",
                   help="Skip writing the B-factor PDB")
    p.add_argument("--no-plot", action="store_true",
                   help="Skip writing the conservation plot")
    return p


def main() -> None:
    args = _build_parser().parse_args()

    if not args.pdb.exists():
        sys.exit(f"[error] PDB not found: {args.pdb}")
    if not args.fasta.exists():
        sys.exit(f"[error] FASTA not found: {args.fasta}")

    print(f"Conservation mapping: {args.pdb.name}  +  {args.fasta.name}")

    try:
        conservation_report(
            pdb_path   = args.pdb,
            fasta_path = args.fasta,
            chain      = args.chain,
            out_dir    = args.out,
            write_pdb  = not args.no_pdb,
            write_plot = not args.no_plot,
        )
    except ValueError as exc:
        sys.exit(f"[error] {exc}")

    print("Done.")


if __name__ == "__main__":
    main()
