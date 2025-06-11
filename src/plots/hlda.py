import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "data/output/FES_HLDA.dat",
    delim_whitespace=True,
    comment="#",
    names=["s_hlda", "F", "dF"],
)

plt.plot(df["s_hlda"], df["F"], marker="o")
plt.xlabel("HLDA coordinate")
plt.ylabel("FES")
plt.tight_layout()
plt.show()
