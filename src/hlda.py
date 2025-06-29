#!/usr/bin/env python3
"""
Compute the HLDA collective variable from two COLVAR files and print:

  • The weight vector W (one weight per descriptor dXX)
  • The first few s_HLDA values for the folded and unfolded sets
  • A convergence plot of HLDA vector as more data is used
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.typing import NDArray
from pathlib import Path

folded_path = Path("data/chignolin/output/COLVAR_PIN")
unfolded_path = Path("data/chignolin/output/COLVAR_FLAT")

def load_descriptors(path: Path, omit: list[str] = []):
    with path.open() as f:
        header_line = next(line for line in f if line.lstrip().startswith("#!"))
        header = header_line.replace("#!", "").strip().split()
    
    df = pd.read_csv(path, comment="#", delim_whitespace=True, names=header, skiprows=1)
    descriptor_cols = [col for col in header if col.startswith("d") and col not in omit]
    
    print(f"{path.name}: using descriptors {descriptor_cols}")
    return df[descriptor_cols].astype(float).to_numpy(), descriptor_cols

omit_descriptors = ["d07", "d08", "d09", "d09", "d15", "d16", "d17", "d18", "d19", "d27"]
A, desc_cols_A = load_descriptors(folded_path, omit_descriptors)
B, desc_cols_B = load_descriptors(unfolded_path, omit_descriptors)

mu_A, mu_B = A.mean(0), B.mean(0)
Sigma_A = np.cov(A, rowvar=False)
Sigma_B = np.cov(B, rowvar=False)
inv_Sigma_A = np.linalg.inv(Sigma_A)
inv_Sigma_B = np.linalg.inv(Sigma_B)
W = (inv_Sigma_A + inv_Sigma_B) @ (mu_A - mu_B)
W /= np.linalg.norm(W)
W_final = W.copy()

print("\nHLDA weights:")
for name, w in zip(desc_cols_A, W):
    print(f"  {name}: {w:+.4f}")

s_folded = A @ W
s_unfolded = B @ W
print("\ns_HLDA (folded):   ", s_folded[:10])
print("\ns_HLDA (unfolded): ", s_unfolded[:10])

Npoints = 20
step_size = min(len(A), len(B)) // Npoints
similarities = []

weights_over_time = []
x = [j * step_size for j in range(1, Npoints + 1)]

for j in range(1, Npoints + 1):
    A_sub = A[:j * step_size]
    B_sub = B[:j * step_size]

    mu_Aj = A_sub.mean(0)
    mu_Bj = B_sub.mean(0)
    Sigma_Aj = np.cov(A_sub, rowvar=False)
    Sigma_Bj = np.cov(B_sub, rowvar=False)

    inv_Sigma_Aj = np.linalg.inv(Sigma_Aj)
    inv_Sigma_Bj = np.linalg.inv(Sigma_Bj)
    Wj = (inv_Sigma_Aj + inv_Sigma_Bj) @ (mu_Aj - mu_Bj)
    Wj /= np.linalg.norm(Wj)

    weights_over_time.append(Wj)
    cos_sim = np.dot(Wj, W_final)
    similarities.append(cos_sim)

weights_over_time = np.array(weights_over_time)
for i, name in enumerate(desc_cols_A):
    plt.plot(x, weights_over_time[:, i], label=name)

plt.legend()
plt.savefig("figures/hlda_convergence_weights.png")
plt.clf()


plt.plot(x, similarities, marker='o')
plt.xlabel("Frames used per class")
plt.ylabel("Cosine similarity")
plt.title("HLDA vector convergence")
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/hlda_convergence_cosine.png")
