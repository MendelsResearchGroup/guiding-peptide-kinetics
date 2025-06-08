import pandas as pd
import matplotlib.pyplot as plt

# Load skipping the comment lines
df = pd.read_csv(
    "data/output/FES_HLDA.dat",
    delim_whitespace=True,
    comment="#",
    names=["s_hlda", "F", "dF"],
)

# Plot the free energy surface
plt.plot(df["s_hlda"], df["F"], marker="o")
plt.xlabel("HLDA coordinate")
plt.ylabel("Free energy (kJ/mol)")
plt.title("Free Energy Surface (FES) along HLDA")
plt.grid(True)
plt.tight_layout()
plt.show()
