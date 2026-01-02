from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from common.utils import _load_colvar_pair, _residues_from_desc


def hlda_from_moments(muA, SA, muB, SB, desc_cols):
    muA = np.asarray(muA, float)
    muB = np.asarray(muB, float)
    SA = np.asarray(SA, float)
    SB = np.asarray(SB, float)

    d = muA - muB

    Sw_inv = np.linalg.inv(SA) + np.linalg.inv(SB)
    Sb = 0.5 * np.outer(d, d)

    eigvals, eigvecs = np.linalg.eig(Sw_inv @ Sb)
    idx = int(np.argmax(np.real(eigvals)))
    lam = float(np.real(eigvals[idx]))
    w = np.real(eigvecs[:, idx])

    Sw = np.linalg.inv(Sw_inv)
    w = w / np.sqrt(float(w.T @ Sw @ w))

    return pd.Series(w, index=desc_cols), lam


def prune(SA, SB, desc_cols, threshold: float):
    """
    Drop highly correlated descriptors separately in folded/unfolded covariances.
    """
    def prune_one(cov):
        keep = list(range(cov.shape[0]))
        while len(keep) > 1:
            sub = cov[np.ix_(keep, keep)]
            std = np.sqrt(np.diag(sub))
            denom = std[:, None] * std[None, :]
            corr = np.where(denom > 0, sub / denom, 0.0)
            
            np.fill_diagonal(corr, 0.0)
            i, j = np.unravel_index(np.argmax(np.abs(corr)), corr.shape)
            if abs(corr[i, j]) < threshold:
                break
            keep.pop(j)
        return set(keep)

    keep_F = prune_one(SA)
    keep_U = prune_one(SB)
    keep_idx = sorted(keep_F & keep_U)
    kept_cols = [desc_cols[i] for i in keep_idx]
    return kept_cols, keep_idx




def compute_lambda_grid(
    base_dir,
    tF_grid: Iterable[float] | None = None,
    tU_grid: Iterable[float] | None = None,
    sample_n: int = 1_000_000,
    n_bins: int = 200,
    rmsd_col: str = "rmsd_ca",
    prune_threshold: float = 0.93,
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
        
        eF, cF, sF, S2F = bin_sufficient_stats(df_F_w, desc_cols, rmsd_col, n_bins)
        eU, cU, sU, S2U = bin_sufficient_stats(df_UF_w, desc_cols, rmsd_col, n_bins)

        cF_cent = centers_from_edges(eF)
        cU_cent = centers_from_edges(eU)

        for tF in tF_grid:
            mask_F = cF_cent <= tF
            muF, covF, nF = aggregate_moments(cF, sF, S2F, mask_F)
            if covF is None or nF < 2:
                continue

            for tU in tU_grid:
                if tU <= tF:
                    continue

                mask_U = cU_cent >= tU
                muU, covU, nU = aggregate_moments(cU, sU, S2U, mask_U)
                if covU is None or nU < 2:
                    continue

                kept_cols, keep_idx = prune(covF, covU, desc_cols, threshold=prune_threshold)

                muF_red = muF[keep_idx]
                covF_red = covF[np.ix_(keep_idx, keep_idx)]
                muU_red = muU[keep_idx]
                covU_red = covU[np.ix_(keep_idx, keep_idx)]

                weights_kept, lam = hlda_from_moments(muF_red, covF_red, muU_red, covU_red, kept_cols)

                full_weights = complete_weights(desc_cols, kept_cols, weights_kept, covF, covU, keep_idx)

                full_weights = {k: round(v, 2) for k, v in full_weights.items()}
                res_weight_sums = [round(v, 2) for v in _aggregate_residue_weights(full_weights)]

                rows.append({
                    "Mutant": protein_dir.name,
                    "tF": float(tF),
                    "tU": float(tU),
                    "lambda": float(lam),
                    "n_desc": len(kept_cols),
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

def _aggregate_residue_weights(weights: dict[str, float], max_res: int = 9) -> list[float]:
    """
    Sum absolute weights per residue index.
    """
    agg = np.zeros(max_res + 1, float)
    for desc, w in weights.items():
        for res_idx in _residues_from_desc(desc):
            if 0 <= res_idx <= max_res:
                agg[res_idx] += abs(float(w))
    return agg.tolist()

__all__ = [
    "hlda_from_moments",
    "prune",
    "complete_weights",
    "bin_sufficient_stats",
    "aggregate_moments",
    "centers_from_edges",
    "moments_from_mask",
    "state_mean",
]
