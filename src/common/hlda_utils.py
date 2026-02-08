from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from common.utils import _load_colvar_pair, _aggregate_residue_weights

from scipy.stats import spearmanr


def hlda_from_moments(muA, SA, muB, SB, desc_cols, ridge=1e-8):
    muA = np.asarray(muA, float)
    muB = np.asarray(muB, float)
    SA = np.asarray(SA, float)
    SB = np.asarray(SB, float)

    d = muA - muB
    Sb = 0.5 * np.outer(d, d)

    Sw_inv = np.linalg.inv(SA) + np.linalg.inv(SB)

    eigvals, eigvecs = np.linalg.eig(Sw_inv @ Sb)
    idx = int(np.argmax(np.real(eigvals)))
    lam = float(np.real(eigvals[idx]))
    w = np.real(eigvecs[:, idx])

    Sw = np.linalg.inv(Sw_inv)
    w = w / np.sqrt(float(w.T @ Sw @ w))

    return pd.Series(w, index=desc_cols), lam


def _spearman_abs_corr(X: np.ndarray) -> np.ndarray:
    if X.size == 0 or X.shape[0] < 2:
        return np.zeros((X.shape[1], X.shape[1]), dtype=float)
    corr, _ = spearmanr(X, axis=0)
    if np.isscalar(corr):
        return np.array([[1.0]], dtype=float)
    corr = np.abs(corr)
    return np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)


def prune(X_F, X_U, desc_cols, threshold: float):
    """
    Drop highly correlated descriptors using Spearman rank correlation computed
    separately on folded/unfolded samples (union of drops).
    """
    corrF = _spearman_abs_corr(X_F)
    corrU = _spearman_abs_corr(X_U)
    corr = np.maximum(corrF, corrU)
    lower = np.tril(corr, k=-1)

    to_drop = [j for j in range(lower.shape[1]) if np.any(lower[:, j] > threshold)]
    keep_idx = [i for i in range(len(desc_cols)) if i not in set(to_drop)]
    kept_cols = [desc_cols[i] for i in keep_idx]
    return kept_cols, keep_idx


def compute_lambda_grid(
    base_dir,
    tF_grid: Iterable[float] | None = None,
    tU_grid: Iterable[float] | None = None,
    sample_n: int = 1000000,
    n_bins: int = 200,
    rmsd_col: str = "rmsd_ca",
    prune_threshold: float = 1,
):
    """
    Compute HLDA eigenvalues for a grid of (tF, tU) RMSD thresholds.

    Returns a DataFrame with one row per mutant / threshold pair, including
    the diagonal variances of folded/unfolded covariances for later Δvar
    calculations.
    """
    tF_grid = np.round(np.linspace(0.18, 0.50, 10), 2) if tF_grid is None else tF_grid
    tU_grid = np.round(np.linspace(0.30, 0.79, 10), 2) if tU_grid is None else tU_grid

    rows = []
    for protein_dir in sorted(Path(base_dir).iterdir()):

        df_F_w, df_UF_w, desc_cols = _load_colvar_pair(protein_dir, sample_n, rmsd_col)
        X_F = df_F_w[desc_cols].to_numpy()
        X_U = df_UF_w[desc_cols].to_numpy()
        rmsd_F = df_F_w[rmsd_col].to_numpy()
        rmsd_U = df_UF_w[rmsd_col].to_numpy()
        eF, cF, sF, S2F = bin_sufficient_stats(df_F_w, desc_cols, rmsd_col, n_bins)
        eU, cU, sU, S2U = bin_sufficient_stats(df_UF_w, desc_cols, rmsd_col, n_bins)

        cF_cent = centers_from_edges(eF)
        cU_cent = centers_from_edges(eU)

        for tF in tF_grid:
            mask_F = cF_cent <= tF
            muF, covF, nF = aggregate_moments(cF, sF, S2F, mask_F)

            for tU in tU_grid:
                if tU <= tF:
                    continue

                mask_U = cU_cent >= tU
                muU, covU, nU = aggregate_moments(cU, sU, S2U, mask_U)

                X_F_thr = X_F[rmsd_F <= tF]
                X_U_thr = X_U[rmsd_U >= tU]
                kept_cols, keep_idx = prune(
                    X_F_thr, X_U_thr, desc_cols, threshold=prune_threshold
                )

                weights_kept, lam = hlda_from_moments(
                    muF[keep_idx],
                    covF[np.ix_(keep_idx, keep_idx)],
                    muU[keep_idx],
                    covU[np.ix_(keep_idx, keep_idx)],
                    kept_cols,
                )

                full_weights = complete_weights(desc_cols, kept_cols, weights_kept, covF, covU, keep_idx)

                full_weights = {k: round(v, 2) for k, v in full_weights.items()}
                res_weight_sums = [round(v, 2) for v in _aggregate_residue_weights(full_weights)]

                rows.append({
                    "Mutant": protein_dir.name,
                    "tF": float(tF),
                    "tU": float(tU),
                    "lambda": float(lam),
                    "nF": nF,
                    "nU": nU,
                    "var_F_diag": np.diag(covF).astype(float).tolist(),
                    "var_U_diag": np.diag(covU).astype(float).tolist(),
                    "weights": full_weights,
                    "res_weights": res_weight_sums,
                })

    return pd.DataFrame(rows)


def load_lambda_grid(
    cache_path: Path | str = Path("../data/hlda_lambda_grid.pkl"),
    force: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """
    Load the cached HLDA (tF, tU) grid or compute it if missing.
    Uses separate Spearman pruning by default.
    """
    cache_path = Path(cache_path)
    if cache_path.exists() and not force:
        return pd.read_pickle(cache_path)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df = compute_lambda_grid(**kwargs)
    df.to_pickle(cache_path)
    return df


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

    mu = subset.mean().to_numpy()
    cov = np.cov(subset.to_numpy(), rowvar=False, ddof=1)
    return mu, cov, int(subset.shape[0])


def complete_weights(desc_cols, kept_cols, weights_kept, covF, covU, keep_idx):
    """
    Extend weights to dropped descriptors by mapping each dropped descriptor
    to the kept descriptor with the strongest correlation.
    """
    full_weights = {name: float(w) for name, w in weights_kept.items()}

    stdF_full = np.sqrt(np.diag(covF))
    stdU_full = np.sqrt(np.diag(covU))
    stdF_kept = stdF_full[keep_idx]
    stdU_kept = stdU_full[keep_idx]
    for j, desc in enumerate(desc_cols):
        if j in keep_idx:
            continue
        denomF = stdF_full[j] * stdF_kept
        denomU = stdU_full[j] * stdU_kept

        corrF = np.where(denomF > 0, covF[j, keep_idx] / denomF, 0.0)
        corrU = np.where(denomU > 0, covU[j, keep_idx] / denomU, 0.0)
        abs_comb = np.maximum(np.abs(corrF), np.abs(corrU))
        best = np.nanargmax(abs_comb)
        mapped_desc = kept_cols[best]
        full_weights[desc] = float(weights_kept.get(mapped_desc, 0.0))

    return full_weights
