"""Zero-shot sequence baseline: ESM-2 masked-marginal log-likelihood ratios.

LLR(variant) = log p(mt | context with the site masked) - log p(wt | same).
Negative values mean the language model finds the substitution unlikely; this is
the standard zero-shot variant-effect score (Meier et al., 2021), used here as
the baseline that structural-confidence metrics have to beat.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .config import TargetConfig
from .curate import load_dataset, load_meta
from .folding import configure_cpu_threads, load_pretrained, resolve_device
from .sequences import AMINO_ACIDS

log = logging.getLogger(__name__)
MAX_RESIDUES = 1022  # ESM-2 positional limit (1024 tokens incl. BOS/EOS)


def _context(cfg: TargetConfig, meta: dict) -> tuple[str, int]:
    """Sequence to score in, and the UniProt number of its first residue."""
    if cfg.llr.context == "construct":
        return meta["construct"]["sequence"], cfg.construct.start
    full = meta["full_sequence"]
    if len(full) <= MAX_RESIDUES:
        return full, 1
    # Longest window that fits, centred on the construct.
    centre = (cfg.construct.start + cfg.construct.end) // 2
    start = min(max(1, centre - MAX_RESIDUES // 2), len(full) - MAX_RESIDUES + 1)
    return full[start - 1:start - 1 + MAX_RESIDUES], start


def compute_llr(cfg: TargetConfig, device: str | None = None, batch_size: int = 8) -> pd.DataFrame:
    import torch
    from transformers import AutoTokenizer, EsmForMaskedLM

    df = load_dataset(cfg)
    variants = df[df["variant_id"] != "WT"]
    context, offset = _context(cfg, load_meta(cfg))
    dev = resolve_device(device or cfg.folding.device)
    if dev.type == "cpu":
        configure_cpu_threads()

    tok = AutoTokenizer.from_pretrained(cfg.llr.model)
    model = load_pretrained(EsmForMaskedLM, cfg.llr.model).eval().to(dev)
    aa_ids = torch.tensor([tok.convert_tokens_to_ids(a) for a in AMINO_ACIDS], device=dev)
    base = tok(context, return_tensors="pt")["input_ids"][0].to(dev)  # <cls> seq <eos>
    if len(base) != len(context) + 2:
        raise RuntimeError("unexpected ESM tokenisation: expected one token per residue")

    for p, wt in zip(variants["position"], variants["wt"]):
        if context[int(p) - offset] != wt:
            raise RuntimeError(f"LLR context has {context[int(p) - offset]} at {p}, dataset says {wt}")
    positions = sorted(int(p) for p in variants["position"].unique())
    log_probs: dict[int, np.ndarray] = {}
    for i in range(0, len(positions), batch_size):
        chunk = positions[i:i + batch_size]
        batch = base.repeat(len(chunk), 1)
        token_idx = torch.tensor([p - offset + 1 for p in chunk], device=dev)
        batch[torch.arange(len(chunk), device=dev), token_idx] = tok.mask_token_id
        with torch.inference_mode():
            logits = model(input_ids=batch).logits
        lp = torch.log_softmax(logits[torch.arange(len(chunk), device=dev), token_idx].float(), dim=-1)
        for p, row in zip(chunk, lp[:, aa_ids].cpu().numpy()):
            log_probs[p] = row
        log.info("LLR: %d/%d positions", min(i + batch_size, len(positions)), len(positions))

    index = {a: k for k, a in enumerate(AMINO_ACIDS)}
    out = pd.DataFrame({
        "variant_id": variants["variant_id"].values,
        "llr": [float(log_probs[int(p)][index[mt]] - log_probs[int(p)][index[wt]])
                for p, wt, mt in zip(variants["position"], variants["wt"], variants["mt"])],
    })
    out["llr_model"] = cfg.llr.model
    out["llr_context"] = f"{cfg.llr.context}:{offset}-{offset + len(context) - 1}"
    out.to_csv(cfg.paths.llr, index=False)
    log.info("wrote %s (%d variants)", cfg.paths.llr, len(out))
    return out
