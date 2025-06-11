import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO

# ---- CONFIG ----
files = ["data/output/COLVAR_PIN", "data/output/COLVAR_FLAT"]
prefix = "d"  # descriptor prefix, like d1, d2, ..., d12

# ---- LOAD & COMBINE ----
dfs = []

for file in files:
    with open(file) as f:
        lines = f.readlines()

    header_line = lines[0]
    if not header_line.startswith("#! FIELDS"):
        raise ValueError(f"Missing '#! FIELDS' line in file: {file}")
    col_names = header_line.strip().split()[2:]
    data_lines = lines[1:]

    df = pd.read_csv(
        StringIO("".join(data_lines)), delim_whitespace=True, names=col_names
    )

    descriptor_cols = [col for col in df.columns if col.startswith(prefix)]
    df_subset = df[descriptor_cols].apply(pd.to_numeric, errors="coerce").dropna()
    dfs.append(df_subset)

# Combine all descriptor data
df_all = pd.concat(dfs, ignore_index=True)

# ---- PLOT ----
corr = df_all.corr()

plt.figure(figsize=(10, 8))
sns.heatmap(corr, annot=True, cmap="coolwarm", vmin=-1, vmax=1, square=True)
plt.title("Correlation Matrix of d1 to d12 (Combined)")
plt.tight_layout()
plt.show()
