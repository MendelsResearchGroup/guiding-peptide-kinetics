from pathlib import Path

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

    names = read_colvar_header_names(f_path)
    desc_cols = [c for c in names if c.startswith("d")]

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

def collect_df(
    all_mfpt,
    mfpt_threshold: float,
    tF: float,
    tU: float,
    skip_short: set[str],
    lambda_df: pd.DataFrame,
    quiet: bool,
) -> pd.DataFrame:
    """
    Assemble MFPT, HLDA eigenvalues, and variance deltas for the given thresholds.

    Parameters
    ----------
    all_mfpt : dict
        Mapping long protein names -> {threshold -> samples}.
    mfpt_threshold : float
        Threshold used for MFPT estimation.
    tF, tU : float
        Folded/unfolded RMSD thresholds used for the HLDA grid.
    skip_short : set[str]
        Short mutant names to skip.
    lambda_df : DataFrame
        Precomputed HLDA grid containing the requested (tF, tU).
    quiet : bool
        If True, suppress per-mutant MFPT output.
    """
    skip_short = set(skip_short)
    required = [p for p in proteins if long_to_short.get(p) not in skip_short]
    missing_proteins = [p for p in required if p not in all_mfpt]
    if missing_proteins:
        raise KeyError(f"Missing MFPT entries for proteins: {', '.join(missing_proteins)}")

    lam_slice = lambda_df[
        np.isclose(lambda_df["tF"], tF) & np.isclose(lambda_df["tU"], tU)
    ].set_index("Mutant")
    if lam_slice.empty:
        raise ValueError(f"lambda_df missing rows for tF={tF}, tU={tU}")
    if "WT" not in lam_slice.index:
        raise KeyError("WT row missing from lambda_df for requested tF/tU")

    wt_vars_F = np.array(lam_slice.loc["WT", "var_F_diag"], float)
    wt_vars_U = np.array(lam_slice.loc["WT", "var_U_diag"], float)

    rows = []
    for long_name in required:
        short = long_to_short.get(long_name)
        medium = short_to_medium.get(short)

        # locate the matching threshold key with tolerance
        thr_key = next(
            (k for k in all_mfpt[long_name].keys() if np.isclose(float(k), mfpt_threshold)),
            None,
        )
        if thr_key is None:
            raise KeyError(f"MFPT threshold {mfpt_threshold} missing for {long_name}")

        mfpt_samples = np.sort(np.array(all_mfpt[long_name][thr_key], float))
        if mfpt_samples.size == 0:
            raise ValueError(f"Empty MFPT samples for {long_name} at threshold {mfpt_threshold}")
        else:
            mfpt, lim = estimateMFPT(mfpt_samples)
            if not quiet:
                print(f"{long_name} ({short}): {(mfpt * 1e-6):.4g} us, extra: {lim:.4g}")

        if short not in lam_slice.index:
            raise KeyError(f"Mutant {short} missing from lambda_df for tF={tF}, tU={tU}")
        lam_row = lam_slice.loc[short]
        varF = np.array(lam_row["var_F_diag"], float)
        varU = np.array(lam_row["var_U_diag"], float)

        abs_dvar_F_sum = float(np.nansum(np.abs(varF - wt_vars_F)))
        abs_dvar_U_sum = float(np.nansum(np.abs(varU - wt_vars_U)))
        rel_dvar_F_sum = float(np.nansum(np.abs(varF - wt_vars_F) / wt_vars_F))
        rel_dvar_U_sum = float(np.nansum(np.abs(varU - wt_vars_U) / wt_vars_U))

        res_weights = lam_row["res_weights"]
        res_weight_cols = {
            f"res_weight_{i}": (res_weights[i] if i < len(res_weights) else np.nan)
            for i in range(10)
        }

        rows.append({
            "long": long_name,
            "medium": medium,
            "short": short,
            "lambda": float(lam_row["lambda"]),
            "tF": tF,
            "tU": tU,
            "mfpt": mfpt,
            "lim": lim,
            "abs_dvar_F": abs_dvar_F_sum,
            "abs_dvar_U": abs_dvar_U_sum,
            "rel_dvar_F": rel_dvar_F_sum,
            "rel_dvar_U": rel_dvar_U_sum,
            "residue_idx": short_to_residue.get(short),
            "property_grp": short_to_property.get(short),
            "Tm": Tm["Tm"].get(short),
            "dTm": (Tm["Tm"].get(short) - Tm["Tm"].get("WT")) if short in Tm.index else np.nan,
            "abs_dTm": abs(Tm["Tm"].get(short) - Tm["Tm"].get("WT")) if short in Tm.index else np.nan,
            "n_desc": int(lam_row["n_desc"]),
            "nF": int(lam_row["nF"]),
            "nU": int(lam_row["nU"]),
            "mean_rmsd_F": float(lam_row["mean_rmsd_F"]) if "mean_rmsd_F" in lam_row else np.nan,
            "mean_rmsd_U": float(lam_row["mean_rmsd_U"]) if "mean_rmsd_U" in lam_row else np.nan,
            **res_weight_cols,
        })

    df = pd.DataFrame(rows)
    df.set_index("short", inplace=True)
    return df
