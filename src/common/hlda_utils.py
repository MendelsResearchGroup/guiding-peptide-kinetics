from __future__ import annotations

import numpy as np
import pandas as pd


def hlda_from_moments(muA, SA, muB, SB, desc_cols, ridge=0.0):
    muA = np.asarray(muA, float)
    muB = np.asarray(muB, float)
    SA = np.asarray(SA, float)
    SB = np.asarray(SB, float)

    p = SA.shape[0]
    if ridge > 0:
        SA = SA + ridge * np.eye(p)
        SB = SB + ridge * np.eye(p)

    d = muA - muB

    uA = np.linalg.solve(SA, d)
    uB = np.linalg.solve(SB, d)
    w = uA + uB

    lam = 0.5 * float(d @ w)

    return pd.Series(w, index=desc_cols), lam


def prune(SA, SB, desc_cols, threshold: float):
    """
    Drop highly correlated descriptors separately in folded/unfolded covariances.
    """
    def prune_one(cov, cols):
        keep = list(range(len(cols)))
        while len(keep) > 1:
            sub_cov = cov[np.ix_(keep, keep)]
            std = np.sqrt(np.diag(sub_cov))
            denom = std[:, None] * std[None, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                corr = np.where(denom > 0, sub_cov / denom, 0.0)
            np.fill_diagonal(corr, 0.0)
            i_max, j_max = np.unravel_index(np.argmax(np.abs(corr)), corr.shape)
            maxcorr = abs(corr[i_max, j_max])
            if maxcorr < threshold:
                break
            keep.pop(j_max)
        return set(keep)

    keep_F = prune_one(SA, desc_cols)
    keep_U = prune_one(SB, desc_cols)
    keep_idx = sorted(list(keep_F & keep_U))
    kept_cols = [desc_cols[i] for i in keep_idx]
    return kept_cols, keep_idx


def bin_sufficient_stats(df: pd.DataFrame, desc_cols: list[str], rmsd_col: str, n_bins: int):
    x = df[rmsd_col].to_numpy()
    X = df[desc_cols].to_numpy()

    edges = np.linspace(x.min(), x.max(), n_bins + 1)
    idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, n_bins - 1)

    p = X.shape[1]
    counts = np.bincount(idx, minlength=n_bins).astype(np.int64)

    sums = np.zeros((n_bins, p), dtype=float)
    for j in range(p):
        sums[:, j] = np.bincount(idx, weights=X[:, j], minlength=n_bins)

    S2 = np.zeros((n_bins, p, p), dtype=float)
    for b in range(n_bins):
        if counts[b] == 0:
            continue
        XB = X[idx == b]
        S2[b] = XB.T @ XB

    return edges, counts, sums, S2


def aggregate_moments(counts, sums, S2, mask: np.ndarray):
    n = counts[mask].sum()
    if n < 2:
        return None, None, int(n)

    s = sums[mask].sum(axis=0)
    S = S2[mask].sum(axis=0)

    mu = s / n
    scatter = S - n * np.outer(mu, mu)
    cov = scatter / (n - 1)

    return mu, cov, int(n)


def centers_from_edges(edges):
    return 0.5 * (edges[:-1] + edges[1:])


def moments_from_mask(df: pd.DataFrame, desc_cols: list[str], mask: pd.Series):
    subset = df.loc[mask, desc_cols]
    if subset.shape[0] < 2:
        return None, None, int(subset.shape[0])
    mu = subset.mean().to_numpy()
    cov = np.cov(subset.to_numpy(), rowvar=False, ddof=1)
    return mu, cov, int(subset.shape[0])


def state_mean(series: pd.Series, mask: pd.Series):
    return float(series.loc[mask].mean()) if mask.any() else np.nan


__all__ = [
    "hlda_from_moments",
    "prune",
    "bin_sufficient_stats",
    "aggregate_moments",
    "centers_from_edges",
    "moments_from_mask",
    "state_mean",
]
