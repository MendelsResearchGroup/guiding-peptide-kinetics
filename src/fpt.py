import pandas as pd
from pathlib import Path

import pandas as pd
from STiMetaD import STiMetaD as STM

base_path = Path("./data/chignolin/output")
max_runs = 50

colnames = ["time", "s_hlda", "metad.bias", "metad.acc", "stop_simulation", "rmsd"]

results = []
warnings = []

for i in range(0, max_runs):
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

    rmsd_threshold = 0.15
    consecutive_count = 2

    condition = df["rmsd"] > rmsd_threshold
    rolling_hits = condition.rolling(consecutive_count).sum() == consecutive_count
    indices = rolling_hits[rolling_hits].index

    if indices.empty:
        warnings.append(
            f"[{index_str}] No {consecutive_count} consecutive frames with rmsd > {rmsd_threshold} found."
        )
        continue

    row = df.loc[indices[0]]
    time = row["time"]
    acc = row["metad.acc"]
    predicted = time * acc
    rmsd = row["rmsd"]

    results.append(
        {"id": i, "time": time, "acc": acc, "predicted": predicted, "rmsd": rmsd}
    )

summary_df = pd.DataFrame(results)
# summary_df.to_csv("folding_summary.csv", index=False)
# print("✅ Done. Results saved to 'folding_summary.csv'.")

if warnings:
    print("\n` Warnings:")
    for w in warnings:
        print(" -", w)

estimator = STM(minSampleSize=5)

if not summary_df.empty:
    samples = summary_df["predicted"].to_numpy()
    mfpt = estimator.estimateMFPT(samples=samples) / 1e6
    rate = estimator.estimateRate(samples=samples) * 1e6
    tstar = estimator.estimateTstar(samples=samples) / 1e6

    print("ST-iMetaD Estimates:")
    print(f"  MFPT  ≈ {mfpt:.6f} μs")
    print(f"  k     ≈ {rate:.6f} 1/μs")
    print(f"  t*    ≈ {tstar:.6f} s")
else:
    print("\n No folding samples collected. Cannot estimate kinetics.")
