"""Minimal PDB handling and structure comparison (NumPy only)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

LDDT_THRESHOLDS = (0.5, 1.0, 2.0, 4.0)
LDDT_RADIUS = 15.0


@dataclass
class CAModel:
    resnum: np.ndarray   # (N,) author residue numbers
    resname: list[str]
    coords: np.ndarray   # (N, 3)
    bfactor: np.ndarray  # (N,)


def read_ca(pdb_text: str, chain: str | None = None) -> CAModel:
    """Cα atoms of the first model (first altloc only, no insertion codes)."""
    resnum, resname, coords, bfac = [], [], [], []
    seen = set()
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        if chain is not None and line[21] != chain:
            continue
        if line[26].strip():  # insertion code
            continue
        n = int(line[22:26])
        if n in seen:          # alternate location
            continue
        seen.add(n)
        resnum.append(n)
        resname.append(line[17:20])
        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        bfac.append(float(line[60:66]) if len(line) >= 66 else 0.0)
    if not coords:
        raise ValueError("no Cα atoms found" + (f" in chain {chain}" if chain else ""))
    return CAModel(np.array(resnum), resname, np.array(coords, dtype=float), np.array(bfac))


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation ``R`` and translation ``t`` minimising ||(mobile @ R.T + t) - target||."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, tc - mc @ r.T


def rmsd(a: np.ndarray, b: np.ndarray, superpose: bool = True) -> float:
    if superpose:
        r, t = kabsch(a, b)
        a = a @ r.T + t
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def lddt_per_residue(model: np.ndarray, reference: np.ndarray,
                     radius: float = LDDT_RADIUS, thresholds=LDDT_THRESHOLDS) -> np.ndarray:
    """Cα-lDDT of ``model`` against ``reference`` (same residue order), per residue.

    Superposition-free: compares all intra-structure distances below ``radius`` in
    the reference, so it isolates *local* rearrangements from rigid-body motion.
    Residues with no neighbours get NaN.
    """
    d_ref = np.linalg.norm(reference[:, None] - reference[None], axis=-1)
    d_mod = np.linalg.norm(model[:, None] - model[None], axis=-1)
    n = len(reference)
    pairs = (d_ref < radius) & ~np.eye(n, dtype=bool)
    diff = np.abs(d_ref - d_mod)
    preserved = np.mean([(diff < t) for t in thresholds], axis=0)
    counts = pairs.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(counts > 0, (preserved * pairs).sum(axis=1) / counts, np.nan)


def lddt(model: np.ndarray, reference: np.ndarray) -> float:
    return float(np.nanmean(lddt_per_residue(model, reference)))


def neighbour_count(coords: np.ndarray, radius: float = 10.0) -> np.ndarray:
    """Cα neighbours within ``radius`` — a simple burial proxy."""
    d = np.linalg.norm(coords[:, None] - coords[None], axis=-1)
    return (d < radius).sum(axis=1) - 1


def transform_pdb(pdb_text: str, r: np.ndarray, t: np.ndarray) -> str:
    """Apply a rigid transform to every ATOM/HETATM record."""
    out = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) @ r.T + t
            line = f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"
        out.append(line)
    return "\n".join(out) + "\n"


def superpose_pdb(pdb_text: str, reference_ca: np.ndarray) -> str:
    """Rigidly move a model onto a reference using all Cα atoms (same residue order)."""
    mobile = read_ca(pdb_text).coords
    if len(mobile) != len(reference_ca):
        raise ValueError("superpose_pdb needs equal-length models")
    r, t = kabsch(mobile, reference_ca)
    return transform_pdb(pdb_text, r, t)
