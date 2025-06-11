import sys
import pandas as pd
import matplotlib.pyplot as plt
import math
from io import StringIO

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
axes = axes.flatten() if n > 1 else [axes]


for i, (file, label) in enumerate(pairs):
    with open(file) as f:
        lines = f.readlines()

    # Extract header from the first line
    header_line = lines[0]
    if header_line.startswith("#! FIELDS"):
        col_names = header_line.strip().split()[2:]  # skip "#!", "FIELDS"
        data_lines = lines[1:]
    else:
        raise ValueError("Missing '#! FIELDS' line in file")

    df = pd.read_csv(
        StringIO("".join(data_lines)), delim_whitespace=True, names=col_names
    )

    # Validate requested label exists
    if label not in df.columns:
        print(f"Label '{label}' not found in file: {file}")
        continue

    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df[label] = pd.to_numeric(df[label], errors="coerce")
    df.dropna(inplace=True)

    x = df["time"].to_numpy()
    y = df[label].to_numpy()

    print(x)
    print(y)
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
