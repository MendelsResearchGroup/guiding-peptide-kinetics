from scipy import optimize, stats
import numpy as np
import pandas as pd
from common.consts import groupByResidue, groupByProperty, long_to_short, short_to_medium, proteins

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


Tm = pd.read_csv('../data/Tm.csv', index_col='Mutant')

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


def collect_df(is_clearer, all_mfpt, th: float):
    avg_change_diff = pd.read_csv(f"../data/average_change_difference{'(clearer state)' if is_clearer else ''}.csv", index_col="Mutant")
    avg_change_folded = pd.read_csv(f"../data/average_change_folded{'(clearer state)' if is_clearer else ''}.csv", index_col="Mutant")
    avg_change_unfolded = pd.read_csv(f"../data/average_change_unfolded{'(clearer state)' if is_clearer else ''}.csv", index_col="Mutant")
    cov_dot_product = pd.read_csv(f"../data/cov_dot_products{'(clearer state)' if is_clearer else ''}.csv", index_col="Mutant")
    variance_differences = pd.read_csv(f"../data/variance_differences{'(clearer state)' if is_clearer else ''}.csv", index_col="Mutant")
    eigenvalue_data = pd.read_csv(f"../data/eigenvalues{'(clearer state)' if is_clearer else ''}.csv", index_col="Mutant")
    
    eigenvalue_data_10_30 = pd.read_csv(f"../data/hlda_summary_10_30.csv", index_col="Mutant")
    eigenvalue_data_12_26 = pd.read_csv(f"../data/hlda_summary_12_26.csv", index_col="Mutant")
    eigenvalue_data_12_30 = pd.read_csv(f"../data/hlda_summary_12_30.csv", index_col="Mutant")
    eigenvalue_data_14_26 = pd.read_csv(f"../data/hlda_summary_14_26.csv", index_col="Mutant")
    eigenvalue_data_14_30 = pd.read_csv(f"../data/hlda_summary_14_30.csv", index_col="Mutant")
    enthalpy = pd.read_csv(f"../data/enthalpy{'(clearer state)' if is_clearer else ''}.csv", index_col="Mutant")
    
    rows = []

    for long_name in proteins:
        short = long_to_short.get(long_name)
        medium = short_to_medium.get(short)

        s = np.sort(np.array(all_mfpt[long_name][th]))
        print(long_name)
        mfpt, lim = estimateMFPT(s)
        print(f"{long_name} ({short}): {(mfpt * 1e-6):.4g} us, extra: {lim:.4g}")

        rows.append({
            "long": long_name,
            medium: medium,
            "short": short,
            "eigenvalue": eigenvalue_data.loc[short, "HLDA_Eigenvalue"],
            # "eigenvalue_10_30": eigenvalue_data_10_30.loc[short, "HLDA_Eigenvalue"],
            # "eigenvalue_12_26": eigenvalue_data_12_26.loc[short, "HLDA_Eigenvalue"],
            # "eigenvalue_12_30": eigenvalue_data_12_30.loc[short, "HLDA_Eigenvalue"],
            # "eigenvalue_14_26": eigenvalue_data_14_26.loc[short, "HLDA_Eigenvalue"],
            # "eigenvalue_14_30": eigenvalue_data_14_30.loc[short, "HLDA_Eigenvalue"],
            "mfpt": mfpt,
            "lim": lim,
            "cossim_F": cov_dot_product['CosSim_Folded'].get(short, None),
            "cossim_U": cov_dot_product['CosSim_Unfolded'].get(short, None),
            "dot_F": cov_dot_product['Dot_Folded'].get(short, None),
            "dot_U": cov_dot_product['Dot_Unfolded'].get(short, None),
            "cos_sim_U": cov_dot_product['CosSim_Unfolded'].get(short, None),
            "abs_dcov_F":  cov_dot_product['abs_dcov_F'].get(short, None),
            "abs_dcov_U":  cov_dot_product['abs_dcov_U'].get(short, None),
            # "avg_change_diff": avg_change_diff['AvgChange_Diff_F-U'].get(short, None),
            "abs_dvar_F": variance_differences['abs_dvar_F'].get(short, None),
            "abs_dvar_U": variance_differences['abs_dvar_U'].get(short, None),   
            "dvar_F(WT-Mut)": variance_differences['dvar_F(WT-Mut)'].get(short, None),
            "dvar_U(WT-Mut)": variance_differences['dvar_U(WT-Mut)'].get(short, None),
            # "avg_change_F": avg_change_folded['AvgChange_Folded'].get(short, None),
            # "avg_change_U": avg_change_unfolded['AvgChange_Unfolded'].get(short, None),
            "enthalpy": enthalpy['dH'].get(short, None),
            "residue_idx": short_to_residue.get(short),
            "property_grp": short_to_property.get(short),
            "Tm": Tm['Tm'].get(short),
        })

    df = pd.DataFrame(rows)
    df.set_index("short", inplace=True)

    return df

def collect_df_best_p(is_clearer, all_mfpt, thresholds):
    suffix = "(clearer state)" if is_clearer else ""
    avg_change_diff   = pd.read_csv(f"../data/average_change_difference{suffix}.csv", index_col="Mutant")
    avg_change_folded = pd.read_csv(f"../data/average_change_folded{suffix}.csv", index_col="Mutant")
    avg_change_unfolded = pd.read_csv(f"../data/average_change_unfolded{suffix}.csv", index_col="Mutant")
    cov_dot_product   = pd.read_csv(f"../data/cov_dot_products{suffix}.csv", index_col="Mutant")
    variance_diffs    = pd.read_csv(f"../data/variance_differences{suffix}.csv", index_col="Mutant")
    eigenvalue_data   = pd.read_csv(f"../data/eigenvalues{suffix}.csv", index_col="Mutant")
    enthalpy          = pd.read_csv(f"../data/enthalpy{suffix}.csv", index_col="Mutant")

    rows = []

    for long_name in proteins:
        short = long_to_short.get(long_name)
        medium = short_to_medium.get(short)

        best_p = -1.0
        best_th = None
        best_mfpt = np.nan
        best_lim = np.nan

        for th in thresholds:
            s = np.sort(np.array(all_mfpt[long_name][th], float))
          
            mfpt, lim = estimateMFPT(s)
            CDF = np.array([i / len(s) for i in range(1, len(s) + 1)])
            fit = optimize.curve_fit(_comulative, s, CDF, p0=(s.mean()))[0]
            fSamples = np.random.exponential(fit[0], size=1000000)
            pval = stats.kstest(s[:lim], fSamples)[1]

            # _, pval, _ = fit_exp_ks(s[:lim])
            
            if np.isfinite(pval) and pval > best_p:
                best_p = pval
                best_th = float(th)
                best_mfpt = mfpt
                best_lim = lim

        rows.append({
            "long": long_name,
            "medium": medium,
            "short": short,
            "best_threshold": best_th,
            "best_p_value": best_p if best_p >= 0 else np.nan,
            "mfpt": best_mfpt,
            "lim": best_lim,
            "eigenvalue": eigenvalue_data.loc[short, "HLDA_Eigenvalue"] if short in eigenvalue_data.index else np.nan,
            "cossim_F": cov_dot_product['CosSim_Folded'].get(short, np.nan),
            "cossim_U": cov_dot_product['CosSim_Unfolded'].get(short, np.nan),
            "dot_F": cov_dot_product['Dot_Folded'].get(short, np.nan),
            "dot_U": cov_dot_product['Dot_Unfolded'].get(short, np.nan),
            "abs_dcov_F": cov_dot_product['abs_dcov_F'].get(short, np.nan),
            "abs_dcov_U": cov_dot_product['abs_dcov_U'].get(short, np.nan),
            "abs_dvar_F": variance_diffs['abs_dvar_F'].get(short, np.nan),
            "abs_dvar_U": variance_diffs['abs_dvar_U'].get(short, np.nan),
            "dvar_F(WT-Mut)": variance_diffs['dvar_F(WT-Mut)'].get(short, np.nan),
            "dvar_U(WT-Mut)": variance_diffs['dvar_U(WT-Mut)'].get(short, np.nan),
            "enthalpy": enthalpy['dH'].get(short, np.nan),
            "residue_idx": short_to_residue.get(short),
            "property_grp": short_to_property.get(short),
            "Tm": Tm['Tm'].get(short, np.nan),
        })

    df = pd.DataFrame(rows).set_index("short")

    # Optional: sort by best_p_value desc, then best_threshold asc
    df = df.sort_values(["best_p_value", "best_threshold"], ascending=[False, True])
    return df
