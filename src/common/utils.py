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