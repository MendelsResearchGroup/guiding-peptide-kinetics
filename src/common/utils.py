from scipy import optimize, stats
import numpy as np
import pandas as pd


def obtainEstimationsDataFrame(samples, minSampleSize):
    samples.sort()

    survival = np.array(
        [(len(samples) - i) / len(samples) for i in range(len(samples))]
    )
    predictions = []
    R2s: list[Unknown] = []
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
    max = data.loc[data.R2 == data.R2.max()]
    limit = max.limit

    return float(1 / max.prediction), int(limit)


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
