"""Benchmark the wild-type prediction against an experimental structure.

Residue correspondence comes from SIFTS (PDBe), so author numbering quirks in the
PDB entry can't silently shift the comparison.
"""

from __future__ import annotations

import json
import logging

import numpy as np
from scipy import stats

from .config import TargetConfig
from .curate import fetch_json, fetch_text, load_meta
from .folding import read_predictions, structure_file
from .structure import kabsch, lddt_per_residue, read_ca

log = logging.getLogger(__name__)

SIFTS_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}"
PDB_URL = "https://files.rcsb.org/download/{pdb}.pdb"
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def author_to_uniprot(sifts: dict, pdb_id: str, accession: str, chain: str) -> dict[int, int]:
    """Map author residue numbers of ``chain`` to UniProt positions."""
    entry = sifts[pdb_id.lower()]["UniProt"].get(accession)
    if entry is None:
        raise ValueError(f"{pdb_id} has no SIFTS mapping to {accession}")
    mapping: dict[int, int] = {}
    for seg in entry["mappings"]:
        if seg.get("chain_id") != chain:
            continue
        a, b = seg["start"].get("author_residue_number"), seg["end"].get("author_residue_number")
        # Unresolved terminal residues have no author number; anchor on whichever end has one.
        # (Every mapped residue is later checked for amino-acid identity, so a wrong offset can't pass.)
        if a is not None:
            offset = seg["unp_start"] - a
        elif b is not None:
            offset = seg["unp_end"] - b
        else:
            continue
        if a is not None and b is not None and seg["unp_end"] - b != offset:
            log.warning("segment %s-%s is not a pure offset; skipping", a, b)
            continue
        for unp in range(seg["unp_start"], seg["unp_end"] + 1):
            mapping[unp - offset] = unp
    if not mapping:
        raise ValueError(f"no SIFTS segments for chain {chain} of {pdb_id}")
    return mapping


def validate_wild_type(cfg: TargetConfig, refresh: bool = False) -> dict | None:
    ref = cfg.reference_structure
    if not ref.pdb_id:
        log.info("no reference structure configured; skipping validation")
        return None
    paths = cfg.paths
    wt_rec = read_predictions(paths.predictions).get("WT")
    if wt_rec is None or not structure_file(paths, wt_rec).exists():
        raise FileNotFoundError("WT prediction missing — run `varfold fold` first")

    pdb_id = ref.pdb_id.upper()
    sifts = fetch_json(SIFTS_URL.format(pdb=pdb_id.lower()), paths.raw / f"sifts_{pdb_id}.json", refresh)
    pdb_text = fetch_text(PDB_URL.format(pdb=pdb_id), paths.raw / f"{pdb_id}.pdb", refresh)
    mapping = author_to_uniprot(sifts, pdb_id, cfg.uniprot, ref.chain)

    exp = read_ca(pdb_text, chain=ref.chain)
    pred = read_ca(structure_file(paths, wt_rec).read_text())
    pred_plddt = np.asarray(wt_rec["plddt"])  # full precision, unlike PDB B-factors
    construct_seq = load_meta(cfg)["construct"]["sequence"]
    start = cfg.construct.start
    pred_by_unp = {start + int(n) - 1: i for i, n in enumerate(pred.resnum)}

    pairs, mismatches = [], 0
    for i, n in enumerate(exp.resnum):
        unp = mapping.get(int(n))
        if unp is None or unp not in pred_by_unp:
            continue
        if THREE_TO_ONE.get(exp.resname[i]) != construct_seq[unp - start]:
            mismatches += 1  # engineered mutations in the crystal construct
            continue
        pairs.append((unp, i, pred_by_unp[unp]))
    if mismatches > 0.1 * (len(pairs) + mismatches):
        raise ValueError(f"{mismatches} residues of {pdb_id} disagree with the sequence — residue mapping is wrong")
    if len(pairs) < 20:
        raise ValueError(f"only {len(pairs)} residues overlap between prediction and {pdb_id}")

    unp_pos = np.array([p[0] for p in pairs])
    x_exp = exp.coords[[p[1] for p in pairs]]
    x_pred = pred.coords[[p[2] for p in pairs]]
    plddt = pred_plddt[[p[2] for p in pairs]]

    r, t = kabsch(x_pred, x_exp)
    dev = np.linalg.norm(x_pred @ r.T + t - x_exp, axis=1)
    per_res_lddt = lddt_per_residue(x_pred, x_exp)
    ok = ~np.isnan(per_res_lddt)
    rho, p_rho = stats.spearmanr(plddt[ok], per_res_lddt[ok])

    result = {
        "pdb_id": pdb_id, "chain": ref.chain, "fingerprint": wt_rec["fingerprint"],
        "n_aligned": len(pairs), "construct_length": len(construct_seq),
        "coverage": round(len(pairs) / len(construct_seq), 3),
        "sequence_mismatches_excluded": mismatches,
        "ca_rmsd": round(float(np.sqrt((dev ** 2).mean())), 3),
        "lddt_ca": round(float(np.nanmean(per_res_lddt)), 4),
        "plddt_vs_lddt_spearman": {"rho": round(float(rho), 3), "p": float(p_rho)},
        "per_residue": [{"position": int(u), "lddt": None if np.isnan(l) else round(float(l), 3),
                         "deviation": round(float(d), 2), "plddt": round(float(b), 2)}
                        for u, l, d, b in zip(unp_pos, per_res_lddt, dev, plddt)],
    }
    paths.validation.write_text(json.dumps(result, indent=2))
    log.info("WT vs %s:%s — %d residues, Cα RMSD %.2f Å, lDDT-Cα %.3f",
             pdb_id, ref.chain, len(pairs), result["ca_rmsd"], result["lddt_ca"])
    return result
