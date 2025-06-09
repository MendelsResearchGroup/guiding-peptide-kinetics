import sys
import pandas as pd
import matplotlib.pyplot as plt
import math

if len(sys.argv) < 3 or (len(sys.argv) - 1) % 2 != 0:
    print("Usage: python plot_property.py file1.xvg label1 file2.xvg label2 ...")
    sys.exit(1)

args: list[str] = sys.argv[1:]
pairs: list[tuple[str, str]] = list(zip(args[::2], args[1::2]))
n: int = len(pairs)

# Auto subplot grid (e.g. 2x2, 3x2, etc.)
cols = min(2, n)
rows = math.ceil(n / cols)

fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4 * rows))
axes = axes.fmlatten() if n > 1 else [axes]

for i, (file, label) in enumerate(pairs):
    df = pd.read_csv(file, sep="\\s+", header=None, names=["time", label])
    x = df["time"].to_numpy()
    y = df[label].to_numpy()

    axes[i].plot(x, y, label=label)
    axes[i].set_xlabel("Time (ps)")
    axes[i].set_ylabel(label)
    axes[i].set_title(f"{label} vs Time")
    axes[i].legend()
    axes[i].grid(True)

# Hide unused subplots if any
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.show()
