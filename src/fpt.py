import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from STiMetaD import STiMetaD 

base_path = Path("./data/chignolin/output")
colnames = ["time", "hlda", "metad.bias", "metad.acc", "stop_simulation", "rmsd"]

rmsd_threshold = 0.1
consecutive_range = range(1, 11)
estimator = STiMetaD(minSampleSize=10)

results_by_count = {k: [] for k in consecutive_range}

for i in range(200):
    index_str = f"{i:03}"
    run_path = base_path / f"run_{index_str}" / f"HLDA_COLVAR_{index_str}"

    if not run_path.exists():
        continue

    try:
        df = pd.read_csv(run_path, sep=r"\s+", comment="#", names=colnames, engine="python")
    except Exception:
        continue

    condition = df["rmsd"] > rmsd_threshold

    for consecutive_count in consecutive_range:
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
            continue

        row = df.iloc[start_index]
        time = row["time"]
        acc = row["metad.acc"]
        predicted = time * acc

        results_by_count[consecutive_count].append(predicted)

mfpt_by_count = []

for k in consecutive_range:
    predicted_time = np.array(results_by_count[k])

    # STiMetaD expects time values in MD steps, which are 2fs in this case
    samples = predicted_time * 500
    mfpt = estimator.estimateMFPT(samples=samples) / 1e6
    mfpt_by_count.append(mfpt)

print(mfpt_by_count)
plt.figure(figsize=(6, 4))
plt.plot(consecutive_range, mfpt_by_count, marker='o')
plt.xlabel("Consecutive Frames Over RMSD Threshold")
plt.ylabel("MFPT (μs)")
plt.title("MFPT vs Consecutive Frame Count")
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/mfpt_vs_consecutive_count.png")
