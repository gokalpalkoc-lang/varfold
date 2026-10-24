"""Statistics for the paper: discrimination, effect sizes, multiple testing, breaking points."""

from __future__ import annotations

import numpy as np
from scipy import stats


def auc(scores: np.ndarray, positives: np.ndarray) -> float:
    """ROC AUC via the Mann–Whitney identity (ties count 1/2)."""
    scores = np.asarray(scores, dtype=float)
    positives = np.asarray(positives, dtype=bool)
    n_pos, n_neg = positives.sum(), (~positives).sum()
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = stats.rankdata(scores)
    return float((ranks[positives].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def roc_curve(scores: np.ndarray, positives: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(FPR, TPR) over all distinct thresholds, from (0,0) to (1,1)."""
    scores = np.asarray(scores, dtype=float)
    positives = np.asarray(positives, dtype=bool)
    order = np.argsort(-scores, kind="mergesort")
    s, y = scores[order], positives[order]
    distinct = np.r_[np.where(np.diff(s))[0], len(s) - 1]
    tps = np.cumsum(y)[distinct]
    fps = (distinct + 1) - tps
    tpr = np.r_[0, tps / max(y.sum(), 1)]
    fpr = np.r_[0, fps / max((~y).sum(), 1)]
    return fpr, tpr


def _stratified_indices(rng: np.random.Generator, positives: np.ndarray, n_boot: int):
    pos_idx, neg_idx = np.where(positives)[0], np.where(~positives)[0]
    for _ in range(n_boot):
        yield np.r_[rng.choice(pos_idx, len(pos_idx)), rng.choice(neg_idx, len(neg_idx))]


def bootstrap_auc(scores, positives, n_boot: int, seed: int) -> tuple[float, float]:
    """Percentile 95% CI, resampling within each class (keeps the class ratio)."""
    scores, positives = np.asarray(scores, float), np.asarray(positives, bool)
    rng = np.random.default_rng(seed)
    vals = [auc(scores[i], positives[i]) for i in _stratified_indices(rng, positives, n_boot)]
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def _delong_components(scores: np.ndarray, positives: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """AUC and its DeLong structural components (Sun & Xu 2014 midrank formulation)."""
    sp, sn = scores[positives], scores[~positives]
    m, n = len(sp), len(sn)
    tz = stats.rankdata(np.r_[sp, sn])
    tx, ty = stats.rankdata(sp), stats.rankdata(sn)
    a = (tz[:m].sum() - m * (m + 1) / 2) / (m * n)
    v10 = (tz[:m] - tx) / n          # per positive: share of negatives it outranks
    v01 = 1 - (tz[m:] - ty) / m      # per negative: share of positives that outrank it
    return float(a), v10, v01


def paired_permutation_p(a: np.ndarray, b: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    """Exact-style paired permutation test of AUC(a) = AUC(b) (Venkatraman & Begg, 1996 style).

    Both scores are rank-transformed, then each variant's pair of ranks is swapped at
    random. Under H0 (the two metrics have the same ROC curve) the swap is exchangeable,
    so the test holds its size even with a small negative class.
    """
    n = len(y)
    ra, rb = stats.rankdata(a) / n, stats.rankdata(b) / n
    m, k = int(y.sum()), int((~y).sum())

    def auc_rows(x):  # AUC for every row of a (P, N) matrix
        ranks = stats.rankdata(x, axis=1)
        return (ranks[:, y].sum(axis=1) - m * (m + 1) / 2) / (m * k)

    observed = abs(auc(ra, y) - auc(rb, y))
    rng = np.random.default_rng(seed)
    hits = 0
    for start in range(0, n_perm, 1000):
        swap = rng.random((min(1000, n_perm - start), n)) < 0.5
        pa, pb = np.where(swap, rb, ra), np.where(swap, ra, rb)
        hits += int((np.abs(auc_rows(pa) - auc_rows(pb)) >= observed - 1e-12).sum())
    return (hits + 1) / (n_perm + 1)


def paired_auc_difference(scores_a, scores_b, positives, n_perm: int = 5000, seed: int = 0) -> dict:
    """AUC(a) − AUC(b) on the same variants.

    Reports the DeLong 95% CI for the difference and two p-values: DeLong's (the
    conventional one) and a paired permutation p (primary). DeLong is liberal with a
    small negative class — 7.0% type-I error at nominal 5% for 305 vs 21 in
    simulation — while the permutation test holds its size.
    """
    a, b, y = np.asarray(scores_a, float), np.asarray(scores_b, float), np.asarray(positives, bool)
    keep = ~(np.isnan(a) | np.isnan(b))  # compare both metrics on exactly the same variants
    a, b, y = a[keep], b[keep], y[keep]
    m, n = int(y.sum()), int((~y).sum())
    if m < 2 or n < 2:
        return {"difference": None, "ci95": None, "p": None, "p_delong": None, "n_pos": m, "n_neg": n}
    auc_a, v10a, v01a = _delong_components(a, y)
    auc_b, v10b, v01b = _delong_components(b, y)
    s10 = np.cov(np.vstack([v10a, v10b]))
    s01 = np.cov(np.vstack([v01a, v01b]))
    s = s10 / m + s01 / n
    var = s[0, 0] + s[1, 1] - 2 * s[0, 1]
    diff = auc_a - auc_b
    if var <= 0:
        p = 1.0 if diff == 0 else 0.0
        half = 0.0
    else:
        se = np.sqrt(var)
        p = float(2 * stats.norm.sf(abs(diff) / se))
        half = 1.959964 * se
    return {"difference": diff, "ci95": [diff - half, diff + half],
            "p": paired_permutation_p(a, b, y, n_perm, seed), "p_delong": p, "n_pos": m, "n_neg": n}


def benjamini_hochberg(pvalues) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def compare_groups(values: np.ndarray, positives: np.ndarray, direction: int, n_boot: int, seed: int) -> dict:
    """Pathogenic (positives) vs benign for one metric.

    ``direction`` orients the metric so larger = more damaging for the AUC; the
    reported medians stay on the metric's own scale.
    """
    raw, y = np.asarray(values, float), np.asarray(positives, bool)
    keep = ~np.isnan(raw)
    raw, y = raw[keep], y[keep]
    if y.sum() < 2 or (~y).sum() < 2:
        return {"n_pos": int(y.sum()), "n_neg": int((~y).sum()), "auc": None}
    s = raw * direction
    a = auc(s, y)
    lo, hi = bootstrap_auc(s, y, n_boot, seed)
    mw = stats.mannwhitneyu(s[y], s[~y], alternative="two-sided")
    return {
        "n_pos": int(y.sum()), "n_neg": int((~y).sum()),
        "auc": a, "auc_ci95": [lo, hi],
        "cliffs_delta": 2 * a - 1,
        "mannwhitney_u": float(mw.statistic), "p": float(mw.pvalue),
        "median_pos": float(np.median(raw[y])), "median_neg": float(np.median(raw[~y])),
    }


def empirical_lower_p(observed: np.ndarray, null: np.ndarray) -> np.ndarray:
    """One-sided P(null <= observed) with the +1 correction (never exactly 0)."""
    null_sorted = np.sort(np.asarray(null, float))
    counts = np.searchsorted(null_sorted, np.asarray(observed, float), side="right")
    return (counts + 1) / (len(null_sorted) + 1)


def fisher_enrichment(in_set: np.ndarray, flagged: np.ndarray) -> dict:
    """2x2 Fisher exact test: is ``flagged`` over-represented within ``in_set``?

    ``odds_ratio`` is None unless finite; ``odds_ratio_kind`` says whether it was
    infinite (no flagged items outside the set) or undefined (nothing flagged at all).
    """
    a = int((in_set & flagged).sum())
    b = int((in_set & ~flagged).sum())
    c = int((~in_set & flagged).sum())
    d = int((~in_set & ~flagged).sum())
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    kind = "finite" if np.isfinite(odds) else ("infinite" if np.isinf(odds) else "undefined")
    return {"table": [[a, b], [c, d]], "odds_ratio": float(odds) if kind == "finite" else None,
            "odds_ratio_kind": kind, "p_one_sided": float(p)}
