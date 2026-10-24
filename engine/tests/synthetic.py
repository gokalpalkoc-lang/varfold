"""Synthetic run directory used by engine and web-backend tests (no model needed)."""

import json

import numpy as np
import pandas as pd

from varfold.config import Construct, TargetConfig
from varfold.folding import fingerprint_dir
from varfold.sequences import AMINO_ACIDS, apply_mutations, Mutation, sequence_hash

THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
         "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
         "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


def make_pdb(seq: str, coords: np.ndarray, bfactors: np.ndarray, chain: str = "A", first: int = 1) -> str:
    lines = []
    for i, (aa, xyz, b) in enumerate(zip(seq, coords, bfactors)):
        lines.append(f"ATOM  {i + 1:5d}  CA  {THREE[aa]} {chain}{first + i:4d}    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{b:6.2f}           C")
    return "\n".join(lines + ["END"]) + "\n"


def random_chain(n: int, rng: np.random.Generator) -> np.ndarray:
    """Compact-ish Cα trace: 3.8 Å steps with a bias back towards the centre."""
    pts = [np.zeros(3)]
    for _ in range(n - 1):
        step = rng.normal(size=3)
        step -= 0.15 * pts[-1] / 3.8
        pts.append(pts[-1] + 3.8 * step / np.linalg.norm(step))
    return np.array(pts)


def build_synthetic_run(tmp_path):
    """A complete fake run directory: dataset, meta, predictions and structures.

    Pathogenic variants get a planted local pLDDT drop and local backbone
    distortion; benign/VUS variants only get small global noise.
    """
    rng = np.random.default_rng(7)
    start, length = 101, 80
    seq = "".join(rng.choice(list(AMINO_ACIDS), length))
    cfg = TargetConfig(name="toy", title="Toy domain", gene="TOY", uniprot="P00000",
                       construct=Construct(start, start + length - 1), output_dir=str(tmp_path))
    cfg.analysis.bootstrap = 200
    paths = cfg.paths
    fp = "test|recycles=4|fp32"
    (paths.structures / fingerprint_dir(fp)).mkdir(parents=True)

    wt_ca = random_chain(length, rng)
    wt_plddt = np.clip(rng.normal(85, 6, length), 40, 98)

    rows = [{"variant_id": "WT", "position": None, "wt": None, "mt": None, "construct_index": None,
             "label": "wild_type", "assertions": "", "review_stars": None, "uniprot_curated": False,
             "clinvar_ids": "", "dbsnp_ids": "", "max_gnomad_af": None, "diseases": "", "in_scan": False,
             "sequence": seq}]
    preds = []

    def add_prediction(vid, sequence, ca, plddt):
        structure = f"{fingerprint_dir(fp)}/{vid}.pdb"
        (paths.structures / structure).write_text(make_pdb(sequence, ca, plddt))
        preds.append({"variant_id": vid, "structure": structure,
                      "sequence_hash": sequence_hash(sequence), "length": len(sequence),
                      "mean_plddt": float(plddt.mean()), "mean_plddt_atoms": float(plddt.mean()),
                      "ptm": 0.8 - 0.002 * float(wt_plddt.mean() - plddt.mean()),
                      "plddt": [round(float(x), 2) for x in plddt], "runtime_s": 1.0, "chunk_size": 64,
                      "device": "cpu (test)", "fingerprint": fp, "model": "test",
                      "num_recycles": 4, "precision": "fp32", "versions": {}, "timestamp": "2026-01-01T00:00:00"})

    add_prediction("WT", seq, wt_ca, wt_plddt)
    used = set()
    labels = ["pathogenic"] * 16 + ["benign"] * 10 + ["vus"] * 6
    for k, label in enumerate(labels):
        while True:
            idx = int(rng.integers(6, length - 6))
            if idx not in used:
                used.add(idx)
                break
        wt_aa = seq[idx]
        mt = next(a for a in rng.permutation(list(AMINO_ACIDS)) if a != wt_aa)
        mut = Mutation(start + idx, wt_aa, mt)
        mseq = apply_mutations(seq, [mut], offset=start)
        ca = wt_ca + rng.normal(0, 0.05, wt_ca.shape)
        plddt = wt_plddt + rng.normal(0, 0.3, length)
        if label == "pathogenic":
            window = slice(idx - 4, idx + 5)
            plddt[window] -= rng.uniform(12, 25)
            ca[window] += rng.normal(0, 1.5, ca[window].shape)
        # Random rigid-body motion: metrics must be invariant to it.
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        q *= np.sign(np.linalg.det(q))
        ca = ca @ q.T + rng.normal(0, 10, 3)
        add_prediction(mut.id, mseq, ca, plddt)
        rows.append({"variant_id": mut.id, "position": mut.position, "wt": wt_aa, "mt": mt,
                     "construct_index": idx + 1, "label": label, "assertions": "P:1" if label == "pathogenic" else "",
                     "review_stars": 1 if k % 2 else None, "uniprot_curated": k % 3 == 0,
                     "clinvar_ids": f"RCV{k:09d}", "dbsnp_ids": "", "max_gnomad_af": None, "diseases": "",
                     "in_scan": False, "sequence": mseq})

    pd.DataFrame(rows).to_csv(paths.dataset, index=False)
    with paths.predictions.open("w") as fh:
        for p in preds:
            fh.write(json.dumps(p) + "\n")
    meta = {"target": {"name": "toy", "title": "Toy domain", "gene": "TOY", "uniprot": "P00000"},
            "retrieved_at": "2026-01-01T00:00:00+00:00", "sources": {}, "filters": {}, "dropped": {},
            "construct": {"start": start, "end": start + length - 1, "length": length, "sequence": seq},
            "full_sequence": "M" * (start - 1) + seq, "annotations": [], "counts": {}}
    paths.dataset_meta.write_text(json.dumps(meta))
    return cfg
