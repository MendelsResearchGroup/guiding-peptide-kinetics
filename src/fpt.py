import pandas as pd
from pathlib import Path

import pandas as pd
from STiMetaD import STiMetaD as STM

base_path = Path("./data/output")
max_runs = 4

colnames = ["time", "s_hlda", "metad.bias", "metad.acc", "stop_simulation", "rmsd"]

results = []
warnings = []

for i in range(max_runs):
    index_str = f"{i:02}"
    run_path = base_path / f"run_0{index_str}" / f"HLDA_COLVAR_0{index_str}"

    if not run_path.exists():
        warnings.append(f"[{index_str}] File not found: {run_path}")
        continue

    try:
        df = pd.read_csv(
            run_path, sep=r"\s+", comment="#", names=colnames, engine="python"
        )
    except Exception as e:
        warnings.append(f"[{index_str}] Error reading {run_path}: {e}")
        continue

    folded_row = df[df["stop_simulation"] == 1.0]
    if folded_row.empty:
        warnings.append(f"[{index_str}] No folding (stop_simulation==1.0) found.")
        continue

    row = folded_row.iloc[0]
    time = row["time"]
    acc = row["metad.acc"]
    predicted = time * acc
    rmsd = row["rmsd"]

    results.append(
        {"id": i, "time": time, "acc": acc, "predicted": predicted, "rmsd": rmsd}
    )

summary_df = pd.DataFrame(results)
summary_df.to_csv("folding_summary.csv", index=False)
print("✅ Done. Results saved to 'folding_summary.csv'.")

if warnings:
    print("\n` Warnings:")
    for w in warnings:
        print(" -", w)

estimator = STM(minSampleSize=5)

if not summary_df.empty:
    samples = summary_df["predicted"].to_numpy()

    mfpt = estimator.estimateMFPT(samples=samples) / 1e6  # convert to µs
    rate = estimator.estimateRate(samples=samples) * 1e6  # convert from 1/ps to 1/µs
    tstar = estimator.estimateTstar(samples=samples) / 1e6

    print("ST-iMetaD Estimates:")
    print(f"  MFPT  ≈ {mfpt:.3f} µs")
    print(f"  k     ≈ {rate:.3f} 1/µs")
    print(f"  t*    ≈ {tstar:.3f} µs")
else:
    print("\n No folding samples collected. Cannot estimate kinetics.")
