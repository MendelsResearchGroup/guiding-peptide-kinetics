#!/usr/bin/env python3
"""
Compute the HLDA collective variable from two COLVAR files and print:

  • The weight vector W (one weight per d1…d12 descriptor)
  • The first few s_HLDA values for the folded and unfolded sets
  • A quick overlap histogram (shown on screen, not saved)

Files:
  COLVAR_PIN   - trajectory trapped in the folded basin
  COLVAR_FLAT  - trajectory trapped in the unfolded basin
"""

import numpy as np
import pandas as pd
from numpy.typing import NDArray

folded_path = "data/output/COLVAR_PIN"
unfolded_path = "data/output/COLVAR_FLAT"

descriptor_cols = [f"d{i}" for i in range(1, 13)]
all_cols = ["time", "rmsd", "e2e", "rg"] + descriptor_cols
ridge = 1e-6


def load_descriptors(path) -> NDArray[np.float64]:
    df = pd.read_csv(
        path, comment="#", sep=r"\s+", names=descriptor_cols, engine="python"
    )
    return df[descriptor_cols].astype(float).to_numpy()


A = load_descriptors(folded_path)
B = load_descriptors(unfolded_path)


mu_A, mu_B = A.mean(0), B.mean(0)
Sigma_A = np.cov(A, rowvar=False)
Sigma_B = np.cov(B, rowvar=False)

# HLDA weight vector (Eq. 6)
inv_Sigma_A = np.linalg.inv(Sigma_A + ridge * np.eye(Sigma_A.shape[0]))
inv_Sigma_B = np.linalg.inv(Sigma_B + ridge * np.eye(Sigma_B.shape[0]))
W = (inv_Sigma_A + inv_Sigma_B) @ (mu_A - mu_B)

W /= np.linalg.norm(W)

# Project every frame → scalar HLDA CV
s_folded = A @ W
s_unfolded = B @ W

print("\nHLDA weights:")
for i, w in enumerate(W, 1):
    print(f"  W_{i} = {w:+.4f}")
print("\ns_HLDA (folded):   ", s_folded[:10])
print("\ns_HLDA (unfolded): ", s_unfolded[:10])


# plt.figure(figsize=(6, 4))
# plt.hist(s_folded, bins=50, density=True, alpha=0.5, label="folded")
# plt.hist(s_unfolded, bins=50, density=True, alpha=0.5, label="unfolded")
# plt.xlabel("s_HLDA"); plt.ylabel("probability density"); plt.legend()
# plt.tight_layout()
# plt.show()
