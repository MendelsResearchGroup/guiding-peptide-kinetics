import numpy as np
import pandas as pd

def obtainEstimationsDataFrame(samples, minSampleSize):
    samples = np.sort(np.asarray(samples, float))
    n = len(samples)
    stop = n // 2
    preds, Ds, limits = [], [], []
    for limit in range(minSampleSize, stop):
        x = samples[:limit]
        tau = x.mean()
        k = 1.0 / tau
        F = 1.0 - np.exp(-x / tau)
        ecdf = np.arange(1, limit + 1) / limit
        D = np.max(np.abs(F - ecdf))
        preds.append(k)
        Ds.append(D)
        limits.append(limit)
    return pd.DataFrame({"prediction": preds, "KS_D": Ds, "limit": limits})

def estimateMFPT(samples, minSampleSize=5):
    df = obtainEstimationsDataFrame(samples, minSampleSize)
    row = df.loc[df.KS_D.idxmin()]
    return float(1.0 / row.prediction), int(row.limit), row.KS_D
