#!/usr/bin/env python3
"""
analyze_ligands.py — Structural interaction analysis for PDB/mmCIF files.

All public functions are importable for use in other projects:

    from analyze_ligands import load_structure, get_ligands
    from analyze_ligands import find_contacts, find_hydrogen_bonds
    from analyze_ligands import find_pi_interactions, find_salt_bridges
    from analyze_ligands import analyze_ligand, to_csv

Selection syntax (PyMOL-compatible subset):
    all, none
    polymer [.protein | .nucleic]
    hetatm, organic, solvent, water
    chain A [+B ...]
    resn LIG [+ATP ...]
    resi 100 | resi 100-200 | resi 100+200
    name CA [+CB ...]
    elem N [+O ...]
    b > 50
    not SELECTION
    SELECTION and SELECTION
    SELECTION or SELECTION
    (SELECTION)
    within N of SELECTION

CLI usage:
    python3 analyze_ligands.py structure.pdb --ligands
    python3 analyze_ligands.py structure.pdb --contacts "resn ATP" "polymer"
    python3 analyze_ligands.py structure.pdb --hbonds "polymer" "resn ATP"
    python3 analyze_ligands.py structure.pdb --pi "polymer" "resn ATP"
    python3 analyze_ligands.py structure.pdb --salt-bridges "polymer" "polymer"
    python3 analyze_ligands.py structure.pdb --analyze "resn ATP" --out report.csv
"""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Optional

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

# ─── Constants ────────────────────────────────────────────────────────────────

WATER_NAMES: frozenset[str] = frozenset({"HOH", "WAT", "H2O", "DOD", "SOL"})

# Common crystallographic excipients — buffers, cryoprotectants, precipitants,
# counter-ions, and detergent fragments seen in PDB HETATM records.
# Used by primary_ligands() to skip junk when no ligand is specified.
COMMON_EXCIPIENTS: frozenset[str] = frozenset({
    # Inorganic ions / salts
    "SO4", "PO4", "HPO", "CL", "NA", "MG", "CA", "ZN", "FE", "MN",
    "CU", "CO", "NI", "K", "IOD", "BR", "F", "CD", "HG", "AU", "PT",
    # Common buffer / precipitant molecules
    "GOL", "EDO", "PEG", "MPD", "DMS", "ACT", "ACY", "EOH", "IPA",
    "EGL", "PGE", "P6G", "1PE", "2PE", "PE4", "PE5", "PE8", "OLC",
    "MES", "EPE", "BIS", "TRS", "HED", "FMT", "ACE", "IMD", "AZI",
    "SCN", "NHE", "NH4", "URE", "GLC", "SUC", "CIT", "TAR", "MLT",
    "BME", "DTT", "TCE", "TBR", "IOH", "MOH", "POL", "DIO", "PDO",
    # Detergents / lipids fragments
    "OLA", "OLB", "PC1", "LMT", "BOG", "DDM", "OG", "PGR",
    # Misc small molecules added during crystallisation
    "NO3", "CO3", "HCO", "FLC", "TLA", "PCA", "CPS",
})

PROTEIN_RESIDUES: frozenset[str] = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    "MSE", "HID", "HIE", "HIP", "CYX", "SEC",
})

NUCLEIC_RESIDUES: frozenset[str] = frozenset({
    "DA", "DT", "DG", "DC", "DU",
    "A", "T", "G", "C", "U",
})

# Aromatic ring atom names for standard amino acids.
# Each residue maps to one or more rings; each ring is a list of atom names.
PROTEIN_AROMATIC_RINGS: dict[str, list[list[str]]] = {
    "PHE": [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
    "TYR": [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
    "TRP": [
        ["CG", "CD1", "NE1", "CE2", "CD2"],           # 5-membered
        ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],   # 6-membered
    ],
    "HIS": [["CG", "ND1", "CD2", "CE1", "NE2"]],
}

# Cationic sidechain atoms (for salt bridges and cation-pi)
CATION_ATOMS: dict[str, list[str]] = {
    "ARG": ["NH1", "NH2", "NE"],
    "LYS": ["NZ"],
    "HIS": ["ND1", "NE2"],
}

# Anionic sidechain atoms
ANION_ATOMS: dict[str, list[str]] = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

# Elements that can act as H-bond donors or acceptors (heavy atoms, no explicit H needed)
HBOND_ELEMENTS: frozenset[str] = frozenset({"N", "O", "S"})

# Default distance cutoffs (Å)
DEFAULT_CONTACT_CUTOFF    = 4.5
DEFAULT_HBOND_CUTOFF      = 3.5
DEFAULT_SALT_CUTOFF       = 4.0
DEFAULT_PI_CUTOFF         = 5.5
DEFAULT_CATION_PI_CUTOFF  = 6.0


# ─── Data Classes ─────────────────────────────────────────────────────────────

@dataclass
class Contact:
    """Generic atom-atom contact."""
    chain1: str
    resn1: str
    resi1: int
    atom1: str
    chain2: str
    resn2: str
    resi2: int
    atom2: str
    distance: float
    interaction_type: str = "contact"

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class HBond:
    """Hydrogen bond detected via heavy-atom geometry (no explicit H required)."""
    donor_chain: str
    donor_resn: str
    donor_resi: int
    donor_atom: str
    acceptor_chain: str
    acceptor_resn: str
    acceptor_resi: int
    acceptor_atom: str
    distance: float    # D···A heavy-atom distance
    angle: float       # pre-D···D···A angle (≈ supplement of D-H···A); 0 if not computable
    interaction_type: str = "hbond"

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class PiInteraction:
    """Pi-stacking (face-to-face / edge-to-face) or cation-pi interaction."""
    chain1: str
    resn1: str
    resi1: int
    ring1_label: str       # comma-separated atom names (or cation atom name)
    chain2: str
    resn2: str
    resi2: int
    ring2_label: str
    center_distance: float
    plane_angle: float     # angle between ring planes (0–90°)
    subtype: str           # "face_to_face" | "edge_to_face" | "intermediate" | "cation_pi"
    interaction_type: str = "pi"

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class SaltBridge:
    """Salt bridge between oppositely charged groups."""
    cation_chain: str
    cation_resn: str
    cation_resi: int
    cation_atom: str
    anion_chain: str
    anion_resn: str
    anion_resi: int
    anion_atom: str
    distance: float
    interaction_type: str = "salt_bridge"

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class LigandReport:
    """All interactions for one ligand residue."""
    ligand_resn: str
    ligand_chain: str
    ligand_resi: int
    contacts: list[Contact]       = field(default_factory=list)
    hbonds: list[HBond]           = field(default_factory=list)
    pi_interactions: list[PiInteraction] = field(default_factory=list)
    salt_bridges: list[SaltBridge]      = field(default_factory=list)

    @property
    def n_contacts(self) -> int: return len(self.contacts)
    @property
    def n_hbonds(self) -> int: return len(self.hbonds)
    @property
    def n_pi(self) -> int: return len(self.pi_interactions)
    @property
    def n_salt_bridges(self) -> int: return len(self.salt_bridges)

    def summary(self) -> str:
        return (
            f"{self.ligand_resn} {self.ligand_chain}:{self.ligand_resi} — "
            f"{self.n_contacts} contacts, {self.n_hbonds} H-bonds, "
            f"{self.n_pi} pi, {self.n_salt_bridges} salt bridges"
        )


# ─── PyMOL-compatible Selection Parser ────────────────────────────────────────

class _TT(Enum):
    KW     = auto()   # keyword (all, not, and, or, within, ...)
    PROP   = auto()   # property (chain, resn, resi, name, elem, b)
    VALUE  = auto()   # bare identifier
    NUMBER = auto()   # numeric literal
    OP     = auto()   # comparison operator
    LPAREN = auto()
    RPAREN = auto()
    PLUS   = auto()
    MINUS  = auto()
    DOT    = auto()

_KEYWORDS   = {"all", "none", "polymer", "hetatm", "organic", "solvent",
               "water", "not", "and", "or", "within", "of", "protein", "nucleic"}
_PROPERTIES = {"chain", "resn", "resi", "name", "elem", "b", "q"}


def _tokenize(s: str) -> list[tuple]:
    tokens: list[tuple] = []
    i = 0
    while i < len(s):
        c = s[i]
        if c.isspace():
            i += 1
        elif c == '(':
            tokens.append((_TT.LPAREN, '(')); i += 1
        elif c == ')':
            tokens.append((_TT.RPAREN, ')')); i += 1
        elif c == '+':
            tokens.append((_TT.PLUS, '+')); i += 1
        elif c == '-':
            tokens.append((_TT.MINUS, '-')); i += 1
        elif c == '.':
            tokens.append((_TT.DOT, '.')); i += 1
        elif s[i:i+2] in ('!=', '<=', '>='):
            tokens.append((_TT.OP, s[i:i+2])); i += 2
        elif c in '<>=':
            tokens.append((_TT.OP, c)); i += 1
        elif c.isdigit():
            j = i
            while j < len(s) and (s[j].isdigit() or s[j] == '.'):
                j += 1
            tokens.append((_TT.NUMBER, float(s[i:j]))); i = j
        elif c.isalpha() or c == '_':
            j = i
            while j < len(s) and (s[j].isalnum() or s[j] in '_'):
                j += 1
            word = s[i:j]
            low = word.lower()
            if low in _KEYWORDS:
                tokens.append((_TT.KW, low))
            elif low in _PROPERTIES:
                tokens.append((_TT.PROP, low))
            else:
                tokens.append((_TT.VALUE, word))
            i = j
        elif c == '*':
            tokens.append((_TT.VALUE, '*')); i += 1
        else:
            i += 1
    return tokens


class SelectionParser:
    """
    Parses PyMOL-compatible selection strings and returns BioPython Atom sets.

    Supported syntax
    ----------------
    all, none
    polymer [.protein | .nucleic]
    hetatm, organic, solvent, water
    chain A [+B ...]
    resn LIG [+ATP ...]
    resi 100 | 100-200 | 100+200
    name CA [+CB ...]
    elem N [+O ...]
    b OP VALUE   (e.g. b > 50)
    not SELECTION
    SELECTION and SELECTION
    SELECTION or SELECTION
    (SELECTION)
    within N of SELECTION
    """

    def __init__(self, structure: Structure):
        self._all: list[Atom] = list(structure.get_atoms())
        self._ns = NeighborSearch(self._all)
        self._tokens: list[tuple] = []
        self._pos: int = 0

    def select_atoms(self, s: str) -> list[Atom]:
        """Return atoms matching the selection string."""
        self._tokens = _tokenize(s.strip())
        self._pos = 0
        return list(self._parse_or())

    def select_residues(self, s: str) -> list[Residue]:
        """Return unique residues containing matching atoms."""
        seen: set[int] = set()
        out: list[Residue] = []
        for a in self.select_atoms(s):
            res = a.get_parent()
            if id(res) not in seen:
                seen.add(id(res))
                out.append(res)
        return out

    # ── internals ────────────────────────────────────────────────────────────

    def _peek(self) -> tuple | None:
        return self._tokens[self._pos] if self._pos < len(self._tokens) else None

    def _consume(self) -> tuple | None:
        tok = self._peek()
        self._pos += 1
        return tok

    def _expect_kw(self, kw: str) -> bool:
        tok = self._peek()
        if tok and tok[0] == _TT.KW and tok[1] == kw:
            self._consume(); return True
        return False

    def _parse_or(self) -> set[Atom]:
        left = self._parse_and()
        while self._peek() == (_TT.KW, 'or'):
            self._consume()
            left |= self._parse_and()
        return left

    def _parse_and(self) -> set[Atom]:
        left = self._parse_not()
        while self._peek() == (_TT.KW, 'and'):
            self._consume()
            left &= self._parse_not()
        return left

    def _parse_not(self) -> set[Atom]:
        if self._peek() == (_TT.KW, 'not'):
            self._consume()
            return set(self._all) - self._parse_not()
        return self._parse_primary()

    def _parse_primary(self) -> set[Atom]:
        tok = self._peek()
        if tok is None:
            return set()
        tt, val = tok

        if tt == _TT.LPAREN:
            self._consume()
            result = self._parse_or()
            if self._peek() and self._peek()[0] == _TT.RPAREN:
                self._consume()
            return result

        if tt == _TT.KW:
            return self._parse_keyword(val)

        if tt == _TT.PROP:
            return self._parse_property(val)

        return set()

    def _parse_keyword(self, kw: str) -> set[Atom]:
        self._consume()

        if kw == 'all':
            return set(self._all)

        if kw == 'none':
            return set()

        if kw == 'polymer':
            poly = PROTEIN_RESIDUES | NUCLEIC_RESIDUES
            if self._peek() and self._peek()[0] == _TT.DOT:
                self._consume()
                sub = self._peek()
                if sub and sub[0] == _TT.KW:
                    self._consume()
                    if sub[1] == 'protein':
                        poly = PROTEIN_RESIDUES
                    elif sub[1] == 'nucleic':
                        poly = NUCLEIC_RESIDUES
            return {a for a in self._all
                    if a.get_parent().get_resname().strip() in poly}

        if kw == 'hetatm':
            return {a for a in self._all if a.get_parent().id[0] not in (' ', '')}

        if kw == 'organic':
            return {a for a in self._all
                    if a.get_parent().id[0] not in (' ', '')
                    and a.get_parent().get_resname().strip() not in WATER_NAMES}

        if kw in ('solvent', 'water'):
            return {a for a in self._all
                    if a.get_parent().get_resname().strip() in WATER_NAMES}

        if kw == 'within':
            n_tok = self._consume()
            dist = float(n_tok[1]) if n_tok else 4.0
            self._expect_kw('of')
            ref_atoms = self._parse_primary()
            result: set[Atom] = set()
            for ref in ref_atoms:
                result.update(self._ns.search(ref.get_vector().get_array(), dist, 'A'))
            return result

        return set()

    def _parse_property(self, prop: str) -> set[Atom]:
        self._consume()

        if prop == 'chain':
            vals = self._parse_str_list()
            return {a for a in self._all if a.get_parent().get_parent().id in vals}

        if prop == 'resn':
            vals = {v.upper() for v in self._parse_str_list()}
            return {a for a in self._all
                    if a.get_parent().get_resname().strip().upper() in vals}

        if prop == 'resi':
            ranges = self._parse_resi_list()
            return {a for a in self._all
                    if _in_ranges(a.get_parent().id[1], ranges)}

        if prop == 'name':
            vals = {v.upper() for v in self._parse_str_list()}
            return {a for a in self._all if a.get_name().strip().upper() in vals}

        if prop == 'elem':
            vals = {v.upper() for v in self._parse_str_list()}
            return {a for a in self._all
                    if a.element and a.element.strip().upper() in vals}

        if prop == 'b':
            op_tok = self._consume()
            val_tok = self._consume()
            if not op_tok or not val_tok:
                return set()
            threshold = float(val_tok[1])
            fn = {'>': float.__gt__, '<': float.__lt__, '=': float.__eq__,
                  '>=': float.__ge__, '<=': float.__le__, '!=': float.__ne__
                  }.get(op_tok[1], lambda a, b: True)
            return {a for a in self._all if fn(a.get_bfactor(), threshold)}

        return set()

    def _parse_str_list(self) -> list[str]:
        vals: list[str] = []

        def _one() -> str | None:
            tok = self._peek()
            if tok is None or tok[0] not in (_TT.VALUE, _TT.NUMBER):
                return None
            self._consume()
            if tok[0] == _TT.NUMBER:
                # Convert float to int string ("2.0" → "2") then check for a
                # trailing VALUE token — handles digit-prefixed residue names
                # like "2HB" or "1PE" that the tokenizer splits into two tokens.
                n = tok[1]
                part = str(int(n)) if n == int(n) else str(n)
                if self._peek() and self._peek()[0] == _TT.VALUE:
                    part += str(self._consume()[1])
                return part
            return str(tok[1])

        val = _one()
        if val is not None:
            vals.append(val)
            while self._peek() and self._peek()[0] == _TT.PLUS:
                self._consume()
                val = _one()
                if val is not None:
                    vals.append(val)
        return vals

    def _parse_resi_list(self) -> list[tuple[int, int]]:
        ranges: list[tuple[int, int]] = []

        def one():
            t = self._peek()
            if t and t[0] == _TT.NUMBER:
                start = int(self._consume()[1])
                if self._peek() and self._peek()[0] == _TT.MINUS:
                    self._consume()
                    end = int(self._consume()[1])
                    ranges.append((start, end))
                else:
                    ranges.append((start, start))
        one()
        while self._peek() and self._peek()[0] == _TT.PLUS:
            self._consume()
            one()
        return ranges


def _in_ranges(n: int, ranges: list[tuple[int, int]]) -> bool:
    return any(lo <= n <= hi for lo, hi in ranges)


# ─── Structure Loading ─────────────────────────────────────────────────────────

def _guess_element(atom_name: str) -> str:
    """Guess element symbol from a PDB atom name when the element column is blank."""
    name = atom_name.strip().lstrip("0123456789")
    if not name:
        return ""
    two = name[:2].upper()
    one = name[0].upper()
    two_char = {"FE", "ZN", "MG", "CA", "MN", "CU", "CO", "NI", "CL", "BR", "SE"}
    return two if two in two_char else one


def load_structure(path: str | Path, structure_id: str = "struct") -> Structure:
    """
    Load a PDB or mmCIF file and return a BioPython Structure object.

    Parameters
    ----------
    path :
        Path to a .pdb, .ent, .cif, or .mmcif file.
    structure_id :
        Label for the structure (default: "struct").
    """
    path = Path(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, str(path))
    # Fill in missing element fields (common for custom/non-standard ligands whose
    # PDB files lack the element column at positions 77-78).
    for atom in structure.get_atoms():
        if not atom.element or not atom.element.strip():
            guessed = _guess_element(atom.get_name())
            if guessed:
                atom.element = guessed
    return structure


# ─── Geometry Utilities ────────────────────────────────────────────────────────

def _coords(atom: Atom) -> np.ndarray:
    return atom.get_vector().get_array()


def _dist(a: Atom, b: Atom) -> float:
    return float(np.linalg.norm(_coords(a) - _coords(b)))


def _ring_center(atoms: list[Atom]) -> np.ndarray:
    return np.mean([_coords(a) for a in atoms], axis=0)


def _ring_normal(atoms: list[Atom]) -> np.ndarray:
    """Unit normal vector to the best-fit plane of the ring atoms."""
    coords = np.array([_coords(a) for a in atoms])
    coords -= coords.mean(axis=0)
    _, _, vh = np.linalg.svd(coords)
    return vh[-1]


def _angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    """Angle in degrees between two vectors (0–180)."""
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom < 1e-10:
        return 0.0
    cos = np.clip(np.dot(v1, v2) / denom, -1.0, 1.0)
    return math.degrees(math.acos(cos))


# ─── Ring Detection ─────────────────────────────────────────────────────────--

@dataclass
class Ring:
    chain: str
    resn: str
    resi: int
    atoms: list[Atom]
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))
    normal: np.ndarray = field(default_factory=lambda: np.zeros(3))
    label: str = ""

    def __post_init__(self):
        if len(self.atoms) >= 3:
            self.center = _ring_center(self.atoms)
            self.normal = _ring_normal(self.atoms)
        self.label = ",".join(a.get_name().strip() for a in self.atoms)


def get_protein_rings(structure: Structure,
                      selection: str = "polymer") -> list[Ring]:
    """
    Return Ring objects for all aromatic residues (PHE/TYR/TRP/HIS) matching selection.

    Parameters
    ----------
    structure :
        BioPython Structure.
    selection :
        PyMOL-like selection to restrict which residues are considered.
    """
    sp = SelectionParser(structure)
    sel_atoms = set(sp.select_atoms(selection))
    rings: list[Ring] = []
    for model in structure:
        for chain in model:
            for res in chain:
                resn = res.get_resname().strip()
                if resn not in PROTEIN_AROMATIC_RINGS:
                    continue
                for ring_names in PROTEIN_AROMATIC_RINGS[resn]:
                    atoms = [res[n] for n in ring_names
                             if n in res and res[n] in sel_atoms]
                    if len(atoms) >= 3:
                        rings.append(Ring(chain=chain.id, resn=resn,
                                         resi=res.id[1], atoms=atoms))
    return rings


def _parse_conect(pdb_path: Path) -> dict[int, set[int]]:
    """Build a serial-number adjacency dict from CONECT records in a PDB file."""
    adj: dict[int, set[int]] = defaultdict(set)
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("CONECT"):
                continue
            parts = line.split()
            try:
                nums = [int(x) for x in parts[1:]]
            except ValueError:
                continue
            src = nums[0]
            for dst in nums[1:]:
                adj[src].add(dst)
                adj[dst].add(src)
    return adj


def _find_rings_in_graph(adj: dict[int, set[int]],
                          min_size: int = 5,
                          max_size: int = 6) -> list[frozenset[int]]:
    """Find unique rings of size min_size..max_size via DFS."""
    found: set[frozenset[int]] = set()

    def dfs(start: int, cur: int, path: list[int], visited: set[int]):
        if len(path) > max_size:
            return
        for nbr in adj[cur]:
            if nbr == start and min_size <= len(path) <= max_size:
                found.add(frozenset(path))
            elif nbr not in visited:
                visited.add(nbr)
                dfs(start, nbr, path + [nbr], visited)
                visited.remove(nbr)

    for node in adj:
        dfs(node, node, [node], {node})
    return list(found)


def get_ligand_rings(structure: Structure,
                     pdb_path: str | Path | None = None) -> list[Ring]:
    """
    Detect aromatic rings in ligand (HETATM) residues using CONECT records.

    Planarity is verified by SVD: rings with a smallest singular value > 0.3 Å
    are discarded as non-planar.

    Parameters
    ----------
    structure :
        BioPython Structure.
    pdb_path :
        Path to the original .pdb file.  Required for CONECT parsing.
        CIF files are not supported (returns empty list).
    """
    if pdb_path is None:
        return []
    pdb_path = Path(pdb_path)
    if pdb_path.suffix.lower() not in ('.pdb', '.ent'):
        return []

    adj = _parse_conect(pdb_path)
    if not adj:
        return []

    # Serial → Atom mapping for HETATM atoms only
    serial_map: dict[int, Atom] = {}
    for atom in structure.get_atoms():
        if atom.get_parent().id[0] not in (' ', ''):
            serial_map[atom.get_serial_number()] = atom

    # Restrict adjacency to HETATM atoms
    hetatm_adj = {
        s: {t for t in nbrs if t in serial_map}
        for s, nbrs in adj.items()
        if s in serial_map
    }

    rings: list[Ring] = []
    for ring_serials in _find_rings_in_graph(hetatm_adj):
        atoms = [serial_map[s] for s in ring_serials if s in serial_map]
        if len(atoms) < 3:
            continue
        # Check planarity via SVD: smallest singular value reflects out-of-plane deviation
        coords = np.array([_coords(a) for a in atoms])
        coords -= coords.mean(axis=0)
        _, sv, _ = np.linalg.svd(coords)
        if sv[-1] > 0.3:
            continue
        res = atoms[0].get_parent()
        rings.append(Ring(chain=res.get_parent().id,
                          resn=res.get_resname().strip(),
                          resi=res.id[1], atoms=atoms))
    return rings


# ─── Core Analysis Functions ───────────────────────────────────────────────────

def get_ligands(structure: Structure,
                exclude_water: bool = True) -> list[Residue]:
    """
    Return all non-polymer HETATM residues.

    Parameters
    ----------
    structure :
        BioPython Structure.
    exclude_water :
        Exclude water molecules (HOH, WAT, etc.). Default True.
    """
    out: list[Residue] = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] == ' ':            # standard polymer residue
                    continue
                resn = res.get_resname().strip()
                if exclude_water and resn in WATER_NAMES:
                    continue
                out.append(res)
    return out


def primary_ligands(structure: Structure,
                    min_atoms: int = 7,
                    exclude_excipients: bool = True) -> list[Residue]:
    """
    Return biologically relevant HETATM residues, largest first.

    Filters out water, common crystallographic excipients, and residues with
    fewer than ``min_atoms`` heavy atoms.  Results are sorted by heavy-atom
    count descending, so index 0 is always the most likely drug/cofactor.

    Parameters
    ----------
    structure :
        BioPython Structure.
    min_atoms :
        Minimum number of heavy atoms required (default 7 — excludes most
        small ions and monoatomic species but keeps e.g. heme, ATP, small
        fragments).
    exclude_excipients :
        Remove residue names in COMMON_EXCIPIENTS (default True).

    Returns
    -------
    List of Residue objects sorted by heavy-atom count descending.
    """
    candidates: list[tuple[int, Residue]] = []
    for res in get_ligands(structure, exclude_water=True):
        resn = res.get_resname().strip()
        if exclude_excipients and resn in COMMON_EXCIPIENTS:
            continue
        n_heavy = sum(1 for a in res if a.element and a.element.strip() != "H")
        if n_heavy < min_atoms:
            continue
        candidates.append((n_heavy, res))
    candidates.sort(key=lambda x: x[0], reverse=True)
    return [r for _, r in candidates]


def find_contacts(
    structure: Structure,
    sel1: str,
    sel2: str = "all",
    cutoff: float = DEFAULT_CONTACT_CUTOFF,
    exclude_same_residue: bool = True,
    exclude_backbone: bool = False,
) -> list[Contact]:
    """
    Find all heavy atom contacts between two selections within a distance cutoff.

    Parameters
    ----------
    structure :
        BioPython Structure.
    sel1, sel2 :
        PyMOL-like selection strings.
    cutoff :
        Maximum distance in Å. Default 4.5 Å.
    exclude_same_residue :
        Skip intra-residue contacts. Default True.
    exclude_backbone :
        Skip backbone atoms (N, CA, C, O). Default False.

    Returns
    -------
    List of Contact objects sorted by distance.
    """
    sp = SelectionParser(structure)
    atoms1 = sp.select_atoms(sel1)
    atoms2_set = set(sp.select_atoms(sel2))
    if not atoms1 or not atoms2_set:
        return []

    backbone = {"N", "CA", "C", "O", "OXT"}
    ns = NeighborSearch(list(atoms2_set))
    seen: set[tuple[int, int]] = set()
    results: list[Contact] = []

    for a1 in atoms1:
        if exclude_backbone and a1.get_name().strip() in backbone:
            continue
        for a2 in ns.search(a1.get_vector().get_array(), cutoff, 'A'):
            if a2 not in atoms2_set:
                continue
            res1, res2 = a1.get_parent(), a2.get_parent()
            if exclude_same_residue and res1 is res2:
                continue
            if exclude_backbone and a2.get_name().strip() in backbone:
                continue
            key = tuple(sorted((id(a1), id(a2))))
            if key in seen:
                continue
            seen.add(key)
            results.append(Contact(
                chain1=res1.get_parent().id, resn1=res1.get_resname().strip(),
                resi1=res1.id[1], atom1=a1.get_name().strip(),
                chain2=res2.get_parent().id, resn2=res2.get_resname().strip(),
                resi2=res2.id[1], atom2=a2.get_name().strip(),
                distance=round(_dist(a1, a2), 3),
            ))

    return sorted(results, key=lambda c: c.distance)


def find_hydrogen_bonds(
    structure: Structure,
    sel1: str = "all",
    sel2: str = "all",
    dist_cutoff: float = DEFAULT_HBOND_CUTOFF,
    angle_cutoff: float = 90.0,
) -> list[HBond]:
    """
    Detect hydrogen bonds between two selections using heavy-atom geometry.

    Since X-ray structures typically lack explicit H positions, bonds are
    identified by:
      1. D···A distance < dist_cutoff (N, O, or S atoms only)
      2. Optionally: pre-D − D − A angle > angle_cutoff (approximates D-H···A)

    Bonds are tested with sel1 as donor and sel2 as acceptor.  Call twice and
    merge to detect bonds in both directions.

    Parameters
    ----------
    structure :
        BioPython Structure.
    sel1 :
        Donor selection.
    sel2 :
        Acceptor selection.
    dist_cutoff :
        Maximum D···A heavy-atom distance (Å). Default 3.5 Å.
    angle_cutoff :
        Minimum pre-D−D−A angle (degrees). Applied only when a pre-donor
        heavy atom is identifiable. Default 90°.

    Returns
    -------
    List of HBond objects sorted by distance.
    """
    sp = SelectionParser(structure)
    donors = [a for a in sp.select_atoms(sel1)
              if a.element and a.element.strip().upper() in HBOND_ELEMENTS
              and not a.get_name().strip().startswith('H')]
    acceptors = [a for a in sp.select_atoms(sel2)
                 if a.element and a.element.strip().upper() in HBOND_ELEMENTS
                 and not a.get_name().strip().startswith('H')]
    if not donors or not acceptors:
        return []

    ns = NeighborSearch(acceptors)
    seen: set[tuple[int, int]] = set()
    results: list[HBond] = []

    for donor in donors:
        for acceptor in ns.search(donor.get_vector().get_array(), dist_cutoff, 'A'):
            if acceptor not in set(acceptors):
                continue
            res_d, res_a = donor.get_parent(), acceptor.get_parent()
            if res_d is res_a:
                continue
            key = tuple(sorted((id(donor), id(acceptor))))
            if key in seen:
                continue
            seen.add(key)

            dist = _dist(donor, acceptor)
            angle = 0.0

            # Approximate D-H···A angle using the heavy atom bonded to D (pre-donor).
            # We look for the nearest heavy atom to the donor within bonding distance.
            pre_candidates = [
                a for a in res_d.get_atoms()
                if a is not donor
                and a.element and a.element.strip().upper() not in ('H', 'D')
                and _dist(donor, a) < 1.9
            ]
            if pre_candidates:
                pre = min(pre_candidates, key=lambda a: _dist(donor, a))
                v_pre = _coords(pre) - _coords(donor)
                v_acc = _coords(acceptor) - _coords(donor)
                angle = _angle_between(v_pre, v_acc)
                if angle < angle_cutoff:
                    continue

            results.append(HBond(
                donor_chain=res_d.get_parent().id,
                donor_resn=res_d.get_resname().strip(),
                donor_resi=res_d.id[1],
                donor_atom=donor.get_name().strip(),
                acceptor_chain=res_a.get_parent().id,
                acceptor_resn=res_a.get_resname().strip(),
                acceptor_resi=res_a.id[1],
                acceptor_atom=acceptor.get_name().strip(),
                distance=round(dist, 3),
                angle=round(angle, 1),
            ))

    return sorted(results, key=lambda h: h.distance)


def find_pi_interactions(
    structure: Structure,
    sel1: str = "polymer",
    sel2: str = "organic",
    pdb_path: str | Path | None = None,
    center_cutoff: float = DEFAULT_PI_CUTOFF,
    cation_pi_cutoff: float = DEFAULT_CATION_PI_CUTOFF,
) -> list[PiInteraction]:
    """
    Detect pi-stacking and cation-pi interactions between two selections.

    Pi-stacking subtypes (ring-center distance < center_cutoff):
      face_to_face   — plane angle  0–30°  (parallel stacking)
      intermediate   — plane angle 30–60°
      edge_to_face   — plane angle 60–90°  (T-shaped / CH-pi)

    Cation-pi: cationic sidechain atom (ARG/LYS/HIS) within cation_pi_cutoff
    of an aromatic ring centroid.

    Ligand aromatic rings are detected from CONECT records when pdb_path is
    supplied (PDB format only; CIF not supported).

    Parameters
    ----------
    structure :
        BioPython Structure.
    sel1, sel2 :
        PyMOL-like selection strings for the two groups.
    pdb_path :
        Path to the original PDB file for ligand ring detection.
    center_cutoff :
        Maximum ring centroid–centroid distance for pi-stacking (Å). Default 5.5 Å.
    cation_pi_cutoff :
        Maximum cation–ring-centroid distance (Å). Default 6.0 Å.

    Returns
    -------
    List of PiInteraction objects sorted by center distance.
    """
    sp = SelectionParser(structure)
    sel1_atoms = set(sp.select_atoms(sel1))
    sel2_atoms = set(sp.select_atoms(sel2))

    lig_rings = get_ligand_rings(structure, pdb_path)

    def rings_for(sel_atoms: set[Atom], sel_str: str) -> list[Ring]:
        protein = get_protein_rings(structure, sel_str)
        lig = [r for r in lig_rings if any(a in sel_atoms for a in r.atoms)]
        return protein + lig

    rings1 = rings_for(sel1_atoms, sel1)
    rings2 = rings_for(sel2_atoms, sel2)

    results: list[PiInteraction] = []
    seen: set[tuple[int, int]] = set()

    # Pi-stacking: ring pairs
    for r1 in rings1:
        for r2 in rings2:
            if r1.chain == r2.chain and r1.resi == r2.resi:
                continue
            key = tuple(sorted((id(r1), id(r2))))
            if key in seen:
                continue
            seen.add(key)

            dist = float(np.linalg.norm(r1.center - r2.center))
            if dist > center_cutoff:
                continue

            angle = _angle_between(r1.normal, r2.normal)
            angle = min(angle, 180.0 - angle)   # fold to 0–90°

            if angle < 30:
                subtype = "face_to_face"
            elif angle < 60:
                subtype = "intermediate"
            else:
                subtype = "edge_to_face"

            results.append(PiInteraction(
                chain1=r1.chain, resn1=r1.resn, resi1=r1.resi, ring1_label=r1.label,
                chain2=r2.chain, resn2=r2.resn, resi2=r2.resi, ring2_label=r2.label,
                center_distance=round(dist, 3),
                plane_angle=round(angle, 1),
                subtype=subtype,
            ))

    # Cation-pi: charged atom near ring centroid
    def _check_cation_pi(rings: list[Ring], cation_sel_atoms: set[Atom]):
        for ring in rings:
            for catom in cation_sel_atoms:
                res_c = catom.get_parent()
                resn_c = res_c.get_resname().strip()
                if resn_c not in CATION_ATOMS:
                    continue
                if catom.get_name().strip() not in CATION_ATOMS[resn_c]:
                    continue
                dist = float(np.linalg.norm(_coords(catom) - ring.center))
                if dist > cation_pi_cutoff:
                    continue
                if res_c.id[1] == ring.resi and res_c.get_parent().id == ring.chain:
                    continue
                results.append(PiInteraction(
                    chain1=res_c.get_parent().id, resn1=resn_c, resi1=res_c.id[1],
                    ring1_label=catom.get_name().strip(),
                    chain2=ring.chain, resn2=ring.resn, resi2=ring.resi,
                    ring2_label=ring.label,
                    center_distance=round(dist, 3),
                    plane_angle=0.0,
                    subtype="cation_pi",
                ))

    _check_cation_pi(rings2, sel1_atoms)
    _check_cation_pi(rings1, sel2_atoms)

    return sorted(results, key=lambda p: p.center_distance)


def find_salt_bridges(
    structure: Structure,
    sel1: str = "polymer",
    sel2: str = "polymer",
    cutoff: float = DEFAULT_SALT_CUTOFF,
) -> list[SaltBridge]:
    """
    Detect salt bridges between oppositely charged groups.

    Cations: ARG (NH1, NH2, NE), LYS (NZ), HIS (ND1, NE2)
    Anions:  ASP (OD1, OD2), GLU (OE1, OE2)

    Works for protein–protein and protein–ligand interactions.
    Ligand residues with standard charged atom names are recognized automatically.

    Parameters
    ----------
    structure :
        BioPython Structure.
    sel1, sel2 :
        PyMOL-like selection strings.  Cation/anion roles are detected
        automatically within each selection.
    cutoff :
        Maximum distance between charged atoms (Å). Default 4.0 Å.

    Returns
    -------
    List of SaltBridge objects sorted by distance.
    """
    sp = SelectionParser(structure)
    atoms1 = set(sp.select_atoms(sel1))
    atoms2 = set(sp.select_atoms(sel2))

    def _charged(atoms: set[Atom], charge_map: dict) -> list[Atom]:
        return [a for a in atoms
                if a.get_parent().get_resname().strip() in charge_map
                and a.get_name().strip() in
                    charge_map[a.get_parent().get_resname().strip()]]

    cations1 = _charged(atoms1, CATION_ATOMS)
    anions1  = _charged(atoms1, ANION_ATOMS)
    cations2 = _charged(atoms2, CATION_ATOMS)
    anions2  = _charged(atoms2, ANION_ATOMS)

    seen: set[tuple[int, int]] = set()
    results: list[SaltBridge] = []

    def _add(cation: Atom, anion: Atom):
        key = tuple(sorted((id(cation), id(anion))))
        if key in seen:
            return
        seen.add(key)
        res_c, res_a = cation.get_parent(), anion.get_parent()
        if res_c is res_a:
            return
        d = _dist(cation, anion)
        if d > cutoff:
            return
        results.append(SaltBridge(
            cation_chain=res_c.get_parent().id,
            cation_resn=res_c.get_resname().strip(),
            cation_resi=res_c.id[1],
            cation_atom=cation.get_name().strip(),
            anion_chain=res_a.get_parent().id,
            anion_resn=res_a.get_resname().strip(),
            anion_resi=res_a.id[1],
            anion_atom=anion.get_name().strip(),
            distance=round(d, 3),
        ))

    for cat in cations1:
        for ani in anions2:
            _add(cat, ani)
    for cat in cations2:
        for ani in anions1:
            _add(cat, ani)

    return sorted(results, key=lambda s: s.distance)


def analyze_ligand(
    structure: Structure,
    ligand_sel: str,
    protein_sel: str = "polymer",
    contact_cutoff: float = DEFAULT_CONTACT_CUTOFF,
    hbond_cutoff: float = DEFAULT_HBOND_CUTOFF,
    salt_bridge_cutoff: float = DEFAULT_SALT_CUTOFF,
    pdb_path: str | Path | None = None,
) -> list[LigandReport]:
    """
    Run all interaction analyses for every ligand residue matching ligand_sel.

    Returns one LigandReport per unique (chain, resn, resi) ligand.

    Parameters
    ----------
    structure :
        BioPython Structure.
    ligand_sel :
        Selection for the ligand(s), e.g. ``"resn ATP"`` or ``"organic"``.
    protein_sel :
        Selection for the binding partner. Default ``"polymer"``.
    contact_cutoff :
        Distance cutoff for non-bonded contacts (Å). Default 4.5 Å.
    hbond_cutoff :
        Distance cutoff for H-bond donor–acceptor distance (Å). Default 3.5 Å.
    salt_bridge_cutoff :
        Distance cutoff for salt bridge charged atoms (Å). Default 4.0 Å.
    pdb_path :
        Path to the original PDB file (enables ligand aromatic ring detection).

    Returns
    -------
    List of LigandReport objects, one per ligand residue.
    """
    sp = SelectionParser(structure)
    lig_residues = sp.select_residues(ligand_sel)
    reports: list[LigandReport] = []

    for res in lig_residues:
        chain = res.get_parent().id
        resn  = res.get_resname().strip()
        resi  = res.id[1]
        this  = f"resn {resn} and chain {chain} and resi {resi}"

        report = LigandReport(ligand_resn=resn, ligand_chain=chain, ligand_resi=resi)

        report.contacts = find_contacts(
            structure, this, protein_sel, cutoff=contact_cutoff)

        # H-bonds in both directions, then deduplicate
        hb_fwd = find_hydrogen_bonds(structure, this, protein_sel, dist_cutoff=hbond_cutoff)
        hb_rev = find_hydrogen_bonds(structure, protein_sel, this, dist_cutoff=hbond_cutoff)
        seen_hb: set[tuple] = set()
        for hb in hb_fwd + hb_rev:
            k = (hb.donor_chain, hb.donor_resi, hb.donor_atom,
                 hb.acceptor_chain, hb.acceptor_resi, hb.acceptor_atom)
            if k not in seen_hb:
                seen_hb.add(k)
                report.hbonds.append(hb)
        report.hbonds.sort(key=lambda h: h.distance)

        report.pi_interactions = find_pi_interactions(
            structure, protein_sel, this, pdb_path=pdb_path)

        report.salt_bridges = find_salt_bridges(
            structure, this, protein_sel, cutoff=salt_bridge_cutoff)

        reports.append(report)

    return reports


# ─── Output Utilities ─────────────────────────────────────────────────────────

def to_csv(interactions: list, path: str | Path) -> None:
    """
    Write a list of interaction dataclass objects to a CSV file.

    Parameters
    ----------
    interactions :
        List of Contact, HBond, PiInteraction, or SaltBridge objects.
    path :
        Output CSV file path.
    """
    if not interactions:
        print("  (no interactions to write)")
        return
    path = Path(path)
    rows = [obj.as_dict() for obj in interactions]
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows)} rows → {path}")


def print_table(interactions: list, max_rows: int = 60) -> None:
    """Pretty-print interactions as a fixed-width table."""
    if not interactions:
        print("  (none)")
        return
    rows = [obj.as_dict() for obj in interactions]
    cols = list(rows[0].keys())
    widths = {c: max(len(c), max(len(str(r[c])) for r in rows)) for c in cols}
    fmt = "  ".join("{:<" + str(widths[c]) + "}" for c in cols)
    print(fmt.format(*cols))
    print("  ".join("-" * widths[c] for c in cols))
    for row in rows[:max_rows]:
        print(fmt.format(*[str(row[c]) for c in cols]))
    if len(rows) > max_rows:
        print(f"  … {len(rows) - max_rows} more rows hidden (use --out to save all)")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Structural interaction analysis for PDB/CIF files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 analyze_ligands.py 1abc.pdb --ligands
  python3 analyze_ligands.py 1abc.pdb --contacts "resn ATP" "polymer"
  python3 analyze_ligands.py 1abc.pdb --hbonds "polymer" "resn ATP"
  python3 analyze_ligands.py 1abc.pdb --pi "polymer" "resn ATP"
  python3 analyze_ligands.py 1abc.pdb --salt-bridges "polymer" "polymer"
  python3 analyze_ligands.py 1abc.pdb --analyze "resn ATP" --out report.csv

PyMOL selection examples:
  "resn ATP"                   residue named ATP
  "chain A"                    chain A only
  "resn ATP and chain B"       ATP in chain B
  "polymer"                    all protein + nucleic residues
  "organic"                    all non-water HETATM residues
  "within 5 of resn ATP"       anything within 5 Å of ATP
  "not solvent"                exclude water
  "resn PHE+TYR+TRP"          phenylalanine, tyrosine, or tryptophan
""",
    )

    parser.add_argument("structure", help="PDB or mmCIF file")

    mode = parser.add_argument_group("analysis mode (pick one)")
    mx = mode.add_mutually_exclusive_group(required=False)
    mx.add_argument("--ligands",      action="store_true",
                    help="List all ligands in the structure")
    mx.add_argument("--contacts",     nargs="+", metavar="SEL",
                    help="Contacts: SEL1 [SEL2]  (default SEL2: all)")
    mx.add_argument("--hbonds",       nargs="+", metavar="SEL",
                    help="H-bonds donor→acceptor: SEL1 [SEL2]")
    mx.add_argument("--pi",           nargs="+", metavar="SEL",
                    help="Pi interactions: SEL1 [SEL2]")
    mx.add_argument("--salt-bridges", nargs="+", metavar="SEL",
                    dest="salt_bridges",
                    help="Salt bridges: SEL1 [SEL2]")
    mx.add_argument("--analyze",      metavar="SEL",
                    help="Full analysis for matching ligand(s). If no mode is "
                         "given, auto-detects primary ligands (non-excipient, "
                         "≥7 heavy atoms) and runs full analysis.")

    params = parser.add_argument_group("parameters")
    params.add_argument("--cutoff",  type=float, default=None,
                        help="Distance cutoff in Å (overrides mode default)")
    params.add_argument("--protein", default="polymer",
                        help="Protein/partner selection for --analyze (default: polymer)")

    out_grp = parser.add_argument_group("output")
    out_grp.add_argument("--out", metavar="FILE",
                          help="Save results to CSV (--analyze saves one file per interaction type)")

    args = parser.parse_args()
    structure = load_structure(args.structure)

    if args.ligands:
        ligs = get_ligands(structure)
        if not ligs:
            print("No ligands found.")
        else:
            print(f"{'Resn':<8} {'Chain':<6} {'Resi':<6} {'Atoms'}")
            print("-" * 32)
            for res in ligs:
                print(f"{res.get_resname().strip():<8} "
                      f"{res.get_parent().id:<6} "
                      f"{res.id[1]:<6} "
                      f"{sum(1 for _ in res.get_atoms())}")

    elif args.contacts:
        sel1 = args.contacts[0]
        sel2 = args.contacts[1] if len(args.contacts) > 1 else "all"
        cutoff = args.cutoff or DEFAULT_CONTACT_CUTOFF
        results = find_contacts(structure, sel1, sel2, cutoff=cutoff)
        print(f"\nContacts  {sel1!r} ↔ {sel2!r}  cutoff={cutoff} Å  →  {len(results)} found\n")
        print_table(results)
        if args.out:
            to_csv(results, args.out)

    elif args.hbonds:
        sel1 = args.hbonds[0]
        sel2 = args.hbonds[1] if len(args.hbonds) > 1 else "all"
        cutoff = args.cutoff or DEFAULT_HBOND_CUTOFF
        results = find_hydrogen_bonds(structure, sel1, sel2, dist_cutoff=cutoff)
        print(f"\nH-bonds  donor={sel1!r}  acceptor={sel2!r}  cutoff={cutoff} Å  →  {len(results)} found\n")
        print_table(results)
        if args.out:
            to_csv(results, args.out)

    elif args.pi:
        sel1 = args.pi[0]
        sel2 = args.pi[1] if len(args.pi) > 1 else "organic"
        results = find_pi_interactions(structure, sel1, sel2, pdb_path=args.structure)
        print(f"\nPi interactions  {sel1!r} ↔ {sel2!r}  →  {len(results)} found\n")
        print_table(results)
        if args.out:
            to_csv(results, args.out)

    elif args.salt_bridges:
        sel1 = args.salt_bridges[0]
        sel2 = args.salt_bridges[1] if len(args.salt_bridges) > 1 else "polymer"
        cutoff = args.cutoff or DEFAULT_SALT_CUTOFF
        results = find_salt_bridges(structure, sel1, sel2, cutoff=cutoff)
        print(f"\nSalt bridges  {sel1!r} ↔ {sel2!r}  cutoff={cutoff} Å  →  {len(results)} found\n")
        print_table(results)
        if args.out:
            to_csv(results, args.out)

    elif args.analyze is not None or not any([
            args.ligands, args.contacts, args.hbonds,
            args.pi, args.salt_bridges]):
        # Explicit selection OR no action given → auto-detect primary ligands
        sel = args.analyze
        if sel is None:
            candidates = primary_ligands(structure)
            if not candidates:
                print("No primary ligands found (all HETATM residues are water, "
                      "excipients, or < 7 heavy atoms).  Use --analyze SEL to "
                      "specify a ligand explicitly.")
                return
            parts = [f"(resn {r.get_resname().strip()} and chain {r.get_parent().id} "
                     f"and resi {r.id[1]})" for r in candidates]
            sel = " or ".join(parts)
            names = ", ".join(
                f"{r.get_resname().strip()} {r.get_parent().id}:{r.id[1]}"
                for r in candidates
            )
            print(f"Auto-detected primary ligand(s): {names}\n")

        reports = analyze_ligand(
            structure, sel,
            protein_sel=args.protein,
            contact_cutoff=args.cutoff or DEFAULT_CONTACT_CUTOFF,
            hbond_cutoff=DEFAULT_HBOND_CUTOFF,
            salt_bridge_cutoff=DEFAULT_SALT_CUTOFF,
            pdb_path=args.structure,
        )
        if not reports:
            print(f"No residues matched selection: {sel!r}")
            return

        for report in reports:
            print(f"\n{'=' * 62}")
            print(f"  {report.summary()}")
            print(f"{'=' * 62}")
            print(f"\nContacts ({report.n_contacts}):")
            print_table(report.contacts)
            print(f"\nH-bonds ({report.n_hbonds}):")
            print_table(report.hbonds)
            print(f"\nPi interactions ({report.n_pi}):")
            print_table(report.pi_interactions)
            print(f"\nSalt bridges ({report.n_salt_bridges}):")
            print_table(report.salt_bridges)

        if args.out:
            stem = Path(args.out).stem
            ext  = Path(args.out).suffix or ".csv"
            for report in reports:
                tag = f"{report.ligand_resn}_{report.ligand_chain}{report.ligand_resi}"
                if report.contacts:
                    to_csv(report.contacts,       f"{stem}_{tag}_contacts{ext}")
                if report.hbonds:
                    to_csv(report.hbonds,         f"{stem}_{tag}_hbonds{ext}")
                if report.pi_interactions:
                    to_csv(report.pi_interactions, f"{stem}_{tag}_pi{ext}")
                if report.salt_bridges:
                    to_csv(report.salt_bridges,   f"{stem}_{tag}_saltbridges{ext}")


if __name__ == "__main__":
    main()
