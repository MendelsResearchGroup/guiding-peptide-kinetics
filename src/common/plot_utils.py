import numpy as np
import pandas as pd


def _default_label(df, label_col):
    if label_col is None:
        return pd.Series(df.index.astype(str), index=df.index)
    return df[label_col].astype(str)


def scatter_with_labels(
    ax,
    df,
    x_col,
    y_col,
    label_col=None,
    color_col="residue_idx",
    colors=None,
    marker_col=None,
    marker_map=None,
    s=70,
    edgecolor="k",
    linewidth=0.6,
    alpha=0.85,
    annotate=True,
):
    """
    Scatter plot with per-row colors and optional labels.

    Parameters
    ----------
    colors : dict
        Mapping from color_col values (int-like) to color strings.
    marker_col : str | None
        Optional column to select markers per row.
    marker_map : dict | None
        Mapping from marker_col values to marker styles.
    """
    if colors is None:
        colors = {}
    if marker_map is None:
        marker_map = {}

    labels = _default_label(df, label_col)
    for i, (_, row) in enumerate(df.iterrows()):
        color_key = row.get(color_col, np.nan)
        try:
            color_key = int(color_key)
        except (TypeError, ValueError):
            color_key = np.nan
        color = colors.get(color_key, "gray")

        marker = "o"
        if marker_col is not None:
            marker = marker_map.get(row.get(marker_col), "o")

        ax.scatter(
            row[x_col],
            row[y_col],
            c=color,
            s=s,
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=alpha,
            marker=marker,
        )
        if annotate:
            ax.text(
                row[x_col],
                row[y_col] + 0.05,
                labels.iloc[i],
                fontsize=8,
                ha="center",
                va="bottom",
            )
