import sys
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

# Check CLI args
if len(sys.argv) < 3 or (len(sys.argv) - 1) % 2 != 0:
    print("Usage: python plot_property_shared.py file1 label1 file2 label2 ...")
    sys.exit(1)

args = sys.argv[1:]
pairs = list(zip(args[::2], args[1::2]))

# Create single plot
fig, ax = plt.subplots(figsize=(8, 5))

for file, label in pairs:
    with open(file) as f:
        lines = f.readlines()

    # Parse header
    header_line = lines[0]
    if not header_line.startswith("#! FIELDS"):
        raise ValueError(f"Missing '#! FIELDS' in file: {file}")
    col_names = header_line.strip().split()[2:]
    data_lines = lines[1:]

    df = pd.read_csv(
        StringIO("".join(data_lines)), delim_whitespace=True, names=col_names
    )

    if label not in df.columns:
        print(f"Label '{label}' not found in {file}")
        continue

    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df[label] = pd.to_numeric(df[label], errors="coerce")
    df.dropna(inplace=True)

    ax.plot(df["time"], df[label], label=label)

ax.set_xlabel("Time (ps)")
ax.set_ylabel("Value")
ax.set_title("Property vs Time")
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()
