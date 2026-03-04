# guiding-peptide-kinetics (Paper Snapshot)

This repository is the fixed snapshot used to reproduce the paper figures and analysis.

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
  - `data/mfpt_all_thresholds-new-ref.pkl`
  - `data/mfpt-pace=25000-new-ref.pkl`
- `hlda_trajectories_*.zip` contains per-mutant trajectory cache data under `data/hlda_trajectories/`

If you need to rebuild archives from an unpacked `data/` tree:

```bash
./scripts/pack_data.sh
```

## MFPT reproducibility options

You can reproduce MFPT-based results in two ways:

1. Generate MFPT samples from FPT simulations using `src/fpt_plumed/` templates (for example through `src/fpt_single_run.sh`).
2. Use the preprocessed MFPT files from `data_core.zip` (recommended for paper reproduction):
   - `data/mfpt_all_thresholds-new-ref.pkl`
   - `data/mfpt-pace=25000-new-ref.pkl`

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
