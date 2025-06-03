import pandas as pd
from pandas import DataFrame

import matplotlib.pyplot as plt

df: DataFrame = pd.read_csv("potential.xvg", sep=r"\s+", header=None, names=["step", "energy"])
df.plot("step")
plt.xlabel("Step")
plt.ylabel("Potential Energy (kJ/mol)")
plt.title("Energy Minimization Potential")
plt.grid(True)
plt.tight_layout()
plt.show()
