import sys
import pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
    print("Usage: python plot_property.py <file> <x_column> <y_column>")
    sys.exit(1)

file, x_col, y_col = sys.argv[1], sys.argv[2], sys.argv[3]

df = pd.read_csv(file, sep='\\s+')

df[x_col] = pd.to_numeric(df[x_col], errors='coerce')
df[y_col] = pd.to_numeric(df[y_col], errors='coerce')
df = df.dropna(subset=[x_col, y_col])

plt.scatter(df[x_col], df[y_col], s=1)
plt.xlabel(x_col)
plt.ylabel(y_col)
plt.tight_layout()
plt.show()
