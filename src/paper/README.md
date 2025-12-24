Paper figures/tables helpers
============================

This folder gathers the plotting/correlation notebooks used for the paper.

- `paper_plots.ipynb` scans MFPT thresholds to pick the best HLDA `(tF, tU)` pair, caches the HLDA grid (`data/hlda_lambda_grid.pkl`), saves the best-threshold summary (`data/hlda_best_thresholds.csv`), exports a per-mutant table for the chosen `(tF, tU)` (`data/hlda_summary_*.csv`), and plots λ correlations versus |ΔTm| and variance deltas.
- `hlda_corr_table.ipynb` (original workflow) shows the fixed-threshold variance and RMSD deltas with λ; kept here for reference alongside the new flow.

Import tip
----------

When running these notebooks from `src/paper/`, ensure the repo root is on `sys.path` so `common.*` imports resolve. In notebooks `__file__` is often undefined, so rely on `Path.cwd()` instead:

```python
import sys
from pathlib import Path
cwd = Path.cwd().resolve()
repo_root = cwd if (cwd / "src").exists() else cwd.parent
sys.path.insert(0, str(repo_root))
```

Workflow
--------

1. Run `paper_plots.ipynb` first; it will compute/cache the HLDA grid if missing and write the summary CSVs under `data/`.
2. Use the exported summary CSVs for downstream statistics/figures, or open `hlda_corr_table.ipynb` if you need the fixed-threshold tables.

Key data products
-----------------

- `data/hlda_lambda_grid.pkl`: cached HLDA grid (λ, variances, mean RMSD per state) across `(tF, tU)`.
- `data/hlda_best_thresholds.csv`: best `(tF, tU)` per MFPT threshold based on |Spearman ρ| (Pearson also reported).
- `data/hlda_summary_tF=*_tU=*.csv`: merged per-mutant table with MFPT, λ, |Δvar| (mean and summed), |ΔTm|, and state RMSDs at the selected `(tF, tU)`.
