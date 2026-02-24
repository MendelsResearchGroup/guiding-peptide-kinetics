import seaborn as sns


def set_paper_style():
    axis_scale = 1.3
    sns.set_theme(
        style="whitegrid",
        context="paper",
        font="DejaVu Sans",
        font_scale=1.0,
        rc={
            # Increase axis-related text by ~30% while keeping paper styling.
            "axes.titlesize": round(11 * axis_scale),
            "axes.labelsize": round(10 * axis_scale),
            "xtick.labelsize": round(9 * axis_scale),
            "ytick.labelsize": round(9 * axis_scale),
            "legend.fontsize": round(9 * axis_scale),
            "legend.title_fontsize": round(9 * axis_scale),
            "figure.titlesize": round(12 * axis_scale),
            # Export defaults for manuscript figures.
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        },
    )
