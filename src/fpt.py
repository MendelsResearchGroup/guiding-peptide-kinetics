import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from STiMetaD import STiMetaD as STM

base_path = Path("./data/chignolin/output")
colnames = ["time", "hlda", "metad.bias", "metad.acc", "stop_simulation", "rmsd"]

results = []
warnings = []

rmsd_threshold = 0.1
consecutive_count = 5
    
for i in range(0, 100):
    index_str = f"{i:03}"
    run_path = base_path / f"run_{index_str}" / f"HLDA_COLVAR_{index_str}"

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

    condition = df["rmsd"] > rmsd_threshold
    streak = 0
    start_index = None

    for idx, val in enumerate(condition):
        if val:
            streak += 1
            if streak == consecutive_count:
                start_index = idx - consecutive_count + 1
                break
        else:
            streak = 0

    if start_index is None:
        warnings.append(
            f"[{index_str}] No {consecutive_count}  consecutive frames with rmsd > {rmsd_threshold} found."
        )
        continue

    row = df.iloc[start_index]

    time = row["time"]
    acc = row["metad.acc"]
    predicted = time * acc
    rmsd = row["rmsd"]

    results.append(
        {"time": time, "acc": acc, "predicted": predicted, "rmsd": rmsd}
    )

summary_df = pd.DataFrame(results)
print(summary_df)
# summary_df.to_csv("folding_summary.csv", index=False)
# print("✅ Done. Results saved to 'folding_summary.csv'.")

if warnings:
    print("\n` Warnings:")
    for w in warnings:
        print(" -", w)

estimator = STM(minSampleSize=10)

if not summary_df.empty:
    samples = summary_df["predicted"].to_numpy() * 500
    
    predicted_time = summary_df["predicted"].to_numpy()
    predicted_time = predicted_time[(predicted_time < 1e6) & (predicted_time > 1000)]
    print(f"mu/sig = {np.mean(predicted_time) / np.std(predicted_time):.4f}")
    plt.hist(predicted_time, bins=30)
    plt.savefig("figures/histogram.png")
    
    mfpt = estimator.estimateMFPT(samples=samples) / 1e6
    rate = estimator.estimateRate(samples=samples) * 1e6
    tstar = estimator.estimateTstar(samples=samples) / 1e6

    print("ST-iMetaD Estimates:")
    print(f"  MFPT  ≈ {mfpt:.6f} μs")
    print(f"  k     ≈ {rate:.6f} 1/μs")
    print(f"  t*    ≈ {tstar:.6f} μs")
else:
    print("\n No folding samples collected. Cannot estimate kinetics.")
