from pathlib import Path
import pandas as pd


def read_colvar_header_names(path: Path) -> list[str]:
    """
    Extract PLUMED column names from the header line.
    """
    with path.open() as f:
        header = next(line for line in f if line.startswith("#!"))
    return header.strip().split()[2:]


def load_colvar(path: Path, usecols: list[str], nrows: int | None = None) -> pd.DataFrame:
    """
    Read a COLVAR file with the given columns and optional row limit.
    """
    return pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        names=read_colvar_header_names(path),
        engine="python",
        usecols=usecols,
        nrows=nrows,
    )


__all__ = ["read_colvar_header_names", "load_colvar"]
