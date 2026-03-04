# guiding-peptide-kinetics (Paper Snapshot)

This repository is the fixed snapshot used to reproduce the paper figures and analysis.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18864706.svg)](https://doi.org/10.5281/zenodo.18864706)

Paper:
- [Guiding Peptide Kinetics via Collective-Variable Tuning of Free-Energy Barriers](https://doi.org/10.48550/arXiv.2602.19936)

## Reproduce the paper flow

1. Create and activate the environment:

```bash
conda env create -f environment.yml
conda activate protein-fes
pip install -e .
```

2. Unpack archived data:

```bash
./scripts/unpack_data.sh
```

3. Run notebooks for paper figures:
- Open notebooks in `src/paper_plots/`
- Execute the required notebooks to regenerate plots

## Data in this repo

Data is stored as split archives in `data_archives/`.

- `data_core.zip` contains shared analysis assets, including preprocessed MFPT files:
  - `data/mfpt_threshold_summaries_ref.pkl`
  - `data/mfpt_samples_pace25000_ref.pkl`
- `hlda_trajectories_*.zip` contains per-mutant trajectory cache data under `data/hlda_trajectories/`

Short description of the MFPT files used in this paper flow:

- `mfpt_samples_pace25000_ref.pkl`: dictionary keyed by mutant, then threshold, containing per-run MFPT samples from the `PACE=25000` setup (typically about 200 runs per mutant/threshold; a few entries are slightly fewer due to missing/failed runs).
- `mfpt_threshold_summaries_ref.pkl`: dictionary keyed by MFPT threshold (`lim`), each value a per-mutant summary DataFrame used by notebooks (for example `mfpt`, `lambda`, `tF`, `tU`, `residue_idx`, `property_grp`, `Tm`, `dTm`, `nF`, `nU`, etc.).
  This summary includes HLDA-derived quantities (for example `lambda`, `tF`, `tU`) through `hlda_lambda_grid.pkl`, which is computed from `data/hlda_trajectories/` via `src/common/hlda_utils.py`.

If you need to rebuild archives from an unpacked `data/` tree:

```bash
./scripts/pack_data.sh
```

## MFPT reproducibility options

You can reproduce MFPT-based results in two ways:

1. Generate MFPT samples from FPT simulations using `src/fpt_plumed/` templates (for example through `src/fpt_single_run.sh`).
2. Use the preprocessed MFPT files from `data_core.zip` (recommended for paper reproduction):
   - `data/mfpt_threshold_summaries_ref.pkl`
   - `data/mfpt_samples_pace25000_ref.pkl`

## How HLDA is run (`hlda_utils`)

HLDA grid generation is implemented in `src/common/hlda_utils.py`.

- `compute_lambda_grid(...)` loads folded/unfolded COLVAR data for each mutant from `data/hlda_trajectories/`, sweeps `(tF, tU)` RMSD thresholds, prunes highly correlated descriptors (Spearman), and computes HLDA weights/eigenvalue (`lambda`) per grid point.
- `load_lambda_grid(...)` is the notebook-facing entrypoint: it loads cached results from `data/hlda_lambda_grid.pkl` if present, otherwise computes and caches them.

Minimal usage pattern (same flow used by paper notebooks):

```python
from pathlib import Path
from common.hlda_utils import load_lambda_grid

data_dir = Path("data")
lambda_grid = load_lambda_grid(
    cache_path=data_dir / "hlda_lambda_grid.pkl",
    base_dir=data_dir / "hlda_trajectories",
    force=False,
)
```

Set `force=True` to recompute the HLDA grid from raw trajectory-derived data.

## Minimal relevant paths

- `src/paper_plots/`: notebooks that generate paper plots
- `src/fpt_plumed/`: PLUMED templates for FPT workflows
- `scripts/unpack_data.sh`: restore `data/` from `*.zip` archives
- `scripts/pack_data.sh`: rebuild split data archives
- `data_archives/`: committed paper snapshot data archives

## Citation

See `CITATION.cff` for software and paper citation metadata.

## License

- Code: MIT (`LICENSE`)
- Paper/manuscript materials: CC BY-NC-ND 4.0 (`LICENSE-paper`)
