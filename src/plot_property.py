# plot_property.py
import sys
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt

if len(sys.argv) < 4 or (len(sys.argv) - 1) % 2 != 0:
    print("Usage: python plot_property.py file1.xvg label1 file2.xvg label2 ...")
    sys.exit(1)

args: list[str] = sys.argv[1:]
pairs: list[tuple[str, str]] = list(zip(args[::2], args[1::2]))

for file, label in pairs:
    df: DataFrame = pd.read_csv(file, sep='\\s+', header=None, names=['time', label])
    
    plt.figure()  # ← create new figure window
    plt.plot(df['time'], df[label], label=label)
    plt.xlabel("Time (ps)")
    plt.ylabel(label.capitalize())
    plt.title(f"{label.capitalize()} vs Time")
    plt.legend()
    plt.tight_layout()

plt.show()
