"""Diagnostic sheet (fiche) of a comparison between a weighting and the default weighting."""
from __future__ import annotations

import numpy as np

__all__ = ["plot_fiche"]


def plot_fiche(res, file=None, title=None):
    """Two-panel diagnostic sheet of a :func:`wecops.compare_weighting` result.

    (a) presence of each labelled space of the default (group of common or proper components
    sharing a label) under the alternative weighting, with the threshold; (b) change (alternative - default) of the four parts of the
    inertia of each block: common, own proper, proper to other blocks, residual.

    Returns the matplotlib figure; saved to ``file`` if given. Requires matplotlib.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    kap = res["kappa"]
    lab = [g["label"] + (f" ({len(g['components'])} comp.)" if len(g["components"]) > 1 else "") for g in res["groups"]]
    D = res["parts_change"]
    blocks = res["block_names"]
    fig, axs = plt.subplots(1, 2, figsize=(12, max(3.5, 0.45 * max(len(lab), len(blocks)) + 1.8)),
                            gridspec_kw=dict(width_ratios=[1, 1.15]))
    ax = axs[0]
    y = np.arange(len(lab))[::-1]
    ax.barh(y, kap, color=["#2a78d6" if k >= res["kappa_min"] else "#eb6834" for k in kap], height=0.6)
    ax.axvline(res["kappa_min"], color="k", ls="--", lw=1)
    ax.set_xlim(0, 1.05)
    ax.set_yticks(y)
    ax.set_yticklabels(lab)
    for yi, k in zip(y, kap):
        ax.text(min(k, 1.0) + 0.01, yi, f"{k:.2f}", va="center", fontsize=9)
    ax.set_xlabel("presence kappa under the alternative weighting")
    verdict = "equivalent" if res["equivalent"] else "interpretation changed"
    ax.set_title(f"(a) {verdict}; motif {'identical' if res['same_motif'] else 'modified'}", loc="left", fontsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax = axs[1]
    cmap = LinearSegmentedColormap.from_list("d", ["#c2185b", "#f7f6f3", "#4e9a2a"])
    ax.imshow(D, cmap=cmap, vmin=-0.3, vmax=0.3, aspect="auto")
    for i in range(D.shape[0]):
        for j in range(4):
            v = 0.0 if abs(D[i, j]) < 0.005 else D[i, j]
            ax.text(j, i, f"{v:+.2f}" if v else "0.00", ha="center", va="center", fontsize=9,
                    fontweight="bold" if abs(v) >= 0.05 else "normal")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["common", "own proper", "proper to\nother blocks", "residual"], fontsize=9)
    ax.set_yticks(range(len(blocks)))
    ax.set_yticklabels(blocks)
    ax.set_title("(b) change of the inertia of each block (alternative - default)", loc="left", fontsize=10)
    fig.suptitle(title or f"Weighting '{res['weighting_name']}' versus default, R = {res['R']}", x=0.01, ha="left")
    fig.tight_layout()
    if file:
        fig.savefig(file, dpi=150, bbox_inches="tight")
    return fig
