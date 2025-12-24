from pathlib import Path
from typing import Iterable

from scipy import optimize, stats
import numpy as np
import pandas as pd

from common.colvar_utils import read_colvar_header_names, load_colvar
from common.consts import (
    groupByResidue,
    groupByProperty,
    long_to_short,
    short_to_medium,
    proteins,
)
from common.hlda_utils import (
    aggregate_moments,
    bin_sufficient_stats,
    centers_from_edges,
    hlda_from_moments,
    prune,
)
import re
import re

def obtainEstimationsDataFrame(samples, minSampleSize):
    samples.sort()

    survival = np.array(
        [(len(samples) - i) / len(samples) for i in range(len(samples))]
    )
    predictions = []
    R2s: list[float] = []
    limits: int = []
    for limit in range(minSampleSize, len(samples)):
        firstSamples = samples[:limit]
        k = -sum(firstSamples * np.log(survival[:limit])) / sum(firstSamples**2)
        predictions.append(k)
        R2s.append(
            1
            - sum((np.log(survival[:limit]) + k * firstSamples) ** 2)
            / sum((np.log(survival[:limit]) - np.log(survival[:limit]).mean()) ** 2)
        )
        limits.append(limit)

    return pd.DataFrame(
        {
            "Tstar": samples[minSampleSize:],
            "prediction": predictions,
            "R2": R2s,
            "limit": limits,
        }
    )

def estimateMFPT(samples, minSampleSize=5):
    """
    Estimates the mean first-passage time.
    """
    data = obtainEstimationsDataFrame(samples=samples, minSampleSize=minSampleSize)
    row = data.loc[data.R2 == data.R2.max()].iloc[0]  # take the first matching row
    limit = int(row.limit)
    mfpt = float(1.0 / row.prediction)
    return mfpt, limit


def _comulative(t, a):
    return 1 - np.exp(-t / a)


def iMetaDMFPT(samples, KStest=False, fitSamples=1000000):
    CDF = np.array([i / len(samples) for i in range(1, len(samples) + 1)])
    samples.sort()
    fit = optimize.curve_fit(_comulative, samples, CDF, p0=(samples.mean()))[0]

    if KStest:
        fSamples = np.random.exponential(fit[0], size=fitSamples)
        pvalue = stats.kstest(samples, fSamples)[1]
        ret = fit[0], pvalue
    else:
        ret = fit[0]

    return ret

short_to_residue = {short: idx for idx, shorts in groupByResidue.items() for short in shorts}
short_to_property = {short: prop for prop, shorts in groupByProperty.items() for short in shorts}

_DATA_DIR = Path(__file__).resolve().parents[2] / "data"
Tm = pd.read_csv(_DATA_DIR / 'Tm.csv', index_col='Mutant')

# wt_Tm = Tm['Tm'].get('WT')

def fit_exp_ks(samples):
    x = np.array(samples, float)
    x.sort()
    n = x.size
    ecdf = np.arange(1, n + 1) / n
    def F(t, tau): return 1 - np.exp(-t / tau)
    tau = optimize.curve_fit(lambda t, tau: F(t, tau), x, ecdf,
                             p0=(x.mean(),), bounds=(0, np.inf))[0][0]
    D, p = stats.kstest(x, "expon", args=(0, tau))
    return tau, p, D


def _mean_abs_diff(a: np.ndarray | None, b: np.ndarray | None) -> float:
    if a is None or b is None:
        return np.nan
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    n = min(a.size, b.size)
    if n == 0:
        return np.nan
    return float(np.nanmean(np.abs(a[:n] - b[:n])))


def _load_colvar_pair(base_dir: Path, sample_n: int, rmsd_col: str):
    f_path = base_dir / "COLVAR_CV_F"
    uf_path = base_dir / "COLVAR_CV_UF"
    if not (f_path.exists() and uf_path.exists()):
        return None, None, None

    names = read_colvar_header_names(f_path)
    desc_cols = [c for c in names if c.startswith("d")]
    if rmsd_col not in names:
        return None, None, None

    usecols = desc_cols + [rmsd_col]
    df_F = load_colvar(f_path, usecols, sample_n)
    df_UF = load_colvar(uf_path, usecols, sample_n)
    return df_F, df_UF, desc_cols


def _residues_from_desc(desc: str) -> list[int]:
    """
    Map a descriptor name like 'd05' or 'd69' to residue indices.
    Assumes each digit corresponds to a residue index (e.g., '05' -> [0,5]).
    """
    m = re.match(r"d(\d+)", desc)
    if not m:
        return []
    digits = m.group(1)
    return [int(ch) for ch in digits]


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



def compute_lambda_grid(
    base_dir: Path = _DATA_DIR / "traj",
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
        if not protein_dir.is_dir():
            continue

        df_F, df_UF, desc_cols = _load_colvar_pair(protein_dir, sample_n, rmsd_col)
        if df_F is None:
            continue

        eF, cF, sF, S2F = bin_sufficient_stats(df_F, desc_cols, rmsd_col, n_bins)
        eU, cU, sU, S2U = bin_sufficient_stats(df_UF, desc_cols, rmsd_col, n_bins)

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

                # simple state RMSD averages (approximate, from bin centers)
                mean_rmsd_F = float((cF_cent[mask_F] * cF[mask_F]).sum() / cF[mask_F].sum()) if cF[mask_F].sum() else np.nan
                mean_rmsd_U = float((cU_cent[mask_U] * cU[mask_U]).sum() / cU[mask_U].sum()) if cU[mask_U].sum() else np.nan

                # Extend weights to all descriptors by mapping dropped descriptors
                # to the kept descriptor with which they are most correlated (folded cov).
                full_weights = {}
                # weights for kept descriptors
                for name, w in weights_kept.items():
                    full_weights[name] = float(w)

                stdF_full = np.sqrt(np.diag(covF))
                stdU_full = np.sqrt(np.diag(covU))
                stdF_kept = stdF_full[keep_idx]
                stdU_kept = stdU_full[keep_idx]
                for j, desc in enumerate(desc_cols):
                    if j in keep_idx:
                        continue
                    denomF = stdF_full[j] * stdF_kept
                    denomU = stdU_full[j] * stdU_kept
                    with np.errstate(divide="ignore", invalid="ignore"):
                        corrF = np.where(denomF > 0, covF[j, keep_idx] / denomF, 0.0)
                        corrU = np.where(denomU > 0, covU[j, keep_idx] / denomU, 0.0)
                    if corrF.size == 0 or (np.all(~np.isfinite(corrF)) and np.all(~np.isfinite(corrU))):
                        continue
                    abs_comb = np.maximum(np.abs(corrF), np.abs(corrU))
                    best = np.nanargmax(abs_comb)
                    mapped_desc = kept_cols[best]
                    full_weights[desc] = float(weights_kept.get(mapped_desc, 0.0))

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
                    "mean_rmsd_F": mean_rmsd_F,
                    "mean_rmsd_U": mean_rmsd_U,
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


def collect_df(
    is_clearer: bool,
    all_mfpt,
    mfpt_threshold: float,
    tF: float | None = 0.25,
    tU: float | None = 0.68,
    skip_short: set[str] | None = None,
    lambda_df: pd.DataFrame | None = None,
    base_dir: Path = Path("../../data/traj"),
    lambda_cache: Path | str = Path("../data/hlda_lambda_grid.pkl"),
    sample_n: int = 1_000_000,
    n_bins: int = 200,
    rmsd_col: str = "rmsd_ca",
    prune_threshold: float = 0.93,
    quiet: bool = False,
) -> pd.DataFrame:
    """
    Assemble MFPT, HLDA eigenvalues, and variance deltas for the given thresholds.

    Parameters
    ----------
    is_clearer : bool
        Ignored (kept for backwards compatibility with notebooks).
    all_mfpt : dict
        Mapping long protein names -> {threshold -> samples}.
    mfpt_threshold : float
        Threshold used for MFPT estimation.
    tF, tU : float, optional
        Folded/unfolded RMSD thresholds used for the HLDA grid. Defaults to
        (0.25, 0.68) for compatibility with older notebooks.
    lambda_df : DataFrame, optional
        Precomputed HLDA grid. If None, it will be loaded or computed.
    """
    tF = 0.25 if tF is None else tF
    tU = 0.68 if tU is None else tU

    lambda_df = lambda_df if lambda_df is not None else load_lambda_grid(
        cache_path=lambda_cache,
        force=False,
        base_dir=base_dir,
        tF_grid=[tF],
        tU_grid=[tU],
        sample_n=sample_n,
        n_bins=n_bins,
        rmsd_col=rmsd_col,
        prune_threshold=prune_threshold,
    )

    lam_slice = lambda_df[
        np.isclose(lambda_df["tF"], tF) & np.isclose(lambda_df["tU"], tU)
    ].set_index("Mutant")

    wt_vars_F = np.array(lam_slice.loc["WT", "var_F_diag"], float) if "WT" in lam_slice.index else None
    wt_vars_U = np.array(lam_slice.loc["WT", "var_U_diag"], float) if "WT" in lam_slice.index else None

    rows = []
    missing = []

    skip_short = set() if skip_short is None else set(skip_short)
    available = {k for k in all_mfpt.keys() if long_to_short.get(k) not in skip_short}
    missing_proteins = [p for p in proteins if p not in available]
    if missing_proteins and not quiet:
        print(f"Warning: missing MFPT entries for proteins: {', '.join(missing_proteins)}")

    for long_name in proteins:
        if long_name not in available:
            continue
        short = long_to_short.get(long_name)
        medium = short_to_medium.get(short)

        # locate the matching threshold key with tolerance
        thr_key = None
        try:
            for k in all_mfpt[long_name].keys():
                if np.isclose(float(k), mfpt_threshold):
                    thr_key = k
                    break
        except Exception:
            thr_key = None

        if thr_key is None:
            missing.append((long_name, mfpt_threshold))
            mfpt_samples = None
        else:
            mfpt_samples = np.sort(np.array(all_mfpt[long_name][thr_key], float))

        if mfpt_samples is None or mfpt_samples.size == 0:
            mfpt = np.nan
            lim = 0
        else:
            mfpt, lim = estimateMFPT(mfpt_samples)
            if not quiet:
                print(f"{long_name} ({short}): {(mfpt * 1e-6):.4g} us, extra: {lim:.4g}")

        lam_row = lam_slice.loc[short] if short in lam_slice.index else None
        varF = np.array(lam_row["var_F_diag"], float) if lam_row is not None else None
        varU = np.array(lam_row["var_U_diag"], float) if lam_row is not None else None

        abs_dvar_F_sum = float(np.nansum(np.abs(varF - wt_vars_F))) if varF is not None and wt_vars_F is not None else np.nan
        abs_dvar_U_sum = float(np.nansum(np.abs(varU - wt_vars_U))) if varU is not None and wt_vars_U is not None else np.nan

        if lam_row is None or "res_weights" not in lam_row:
            raise ValueError("HLDA residue weights missing; rebuild hlda_lambda_grid.pkl with weights.")
        res_weights = lam_row["res_weights"]
        res_weight_cols = {
            f"res_weight_{i}": (res_weights[i] if i < len(res_weights) else np.nan)
            for i in range(10)
        }

        rows.append({
            "long": long_name,
            "medium": medium,
            "short": short,
            "lambda": float(lam_row["lambda"]) if lam_row is not None else np.nan,
            "tF": tF,
            "tU": tU,
            "mfpt": mfpt,
            "lim": lim,
            "abs_dvar_F": _mean_abs_diff(varF, wt_vars_F),
            "abs_dvar_U": _mean_abs_diff(varU, wt_vars_U),
            "abs_dvar_mean": np.nanmean([_mean_abs_diff(varF, wt_vars_F), _mean_abs_diff(varU, wt_vars_U)]),
            "abs_dvar_F_sum": abs_dvar_F_sum,
            "abs_dvar_U_sum": abs_dvar_U_sum,
            "abs_dvar_sum_mean": np.nanmean([abs_dvar_F_sum, abs_dvar_U_sum]),
            "residue_idx": short_to_residue.get(short),
            "property_grp": short_to_property.get(short),
            "Tm": Tm["Tm"].get(short),
            "dTm": (Tm["Tm"].get(short) - Tm["Tm"].get("WT")) if short in Tm.index else np.nan,
            "abs_dTm": abs(Tm["Tm"].get(short) - Tm["Tm"].get("WT")) if short in Tm.index else np.nan,
            "n_desc": int(lam_row["n_desc"]) if lam_row is not None else np.nan,
            "nF": int(lam_row["nF"]) if lam_row is not None else np.nan,
            "nU": int(lam_row["nU"]) if lam_row is not None else np.nan,
            "mean_rmsd_F": float(lam_row["mean_rmsd_F"]) if lam_row is not None and "mean_rmsd_F" in lam_row else np.nan,
            "mean_rmsd_U": float(lam_row["mean_rmsd_U"]) if lam_row is not None and "mean_rmsd_U" in lam_row else np.nan,
            **res_weight_cols,
        })

    df = pd.DataFrame(rows)
    df.set_index("short", inplace=True)
    if missing:
        details = ", ".join(f"{name}@{thr}" for name, thr in missing)
        raise ValueError(f"MFPT samples missing for: {details}")
    return df
