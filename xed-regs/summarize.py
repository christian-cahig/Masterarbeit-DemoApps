"""
A script summarizing the results of varying the supply regularization strengths

Masterarbeit-DemApp
Copyright (C) 2022 - present, Christian Cahig
https://github.com/christian-cahig/Masterarbeit-DemApp

This file is part of Masterarbeit-DemApp, and is covered by the CC-BY-4.0 License.
See `LICENSE` file in the source directory for details.
"""

from argparse import ArgumentParser

import numpy as np
from numpy.typing import NDArray
from scipy.io import loadmat
import matplotlib.pyplot as plt

__author__ = "Christian Cahig"
__version__ = "0.1.0"
__all__ = []

SYSTMS = ["ACT500", "RTE1888", "POL3375"]

def get_argparser() -> ArgumentParser:
    """
    Returns the argument parser
    """
    ap = ArgumentParser(description="A script for summarizing the results")
    ap.add_argument(
        "--system", type=str, required=False, default="ACT500", choices=SYSTMS,
        help="System whose results are to be summarized",
    )
    ap.add_argument(
        "--slack_scale", type=int, required=False, default=None,
        help="Scales the distributed slack quantities by 10^(-SLACK_SCALE)",
    )
    ap.add_argument(
        "--marker_face_color", type=str, required=False, default="blue",
        help="Matplotlib-supported named color for the marker faces",
    )
    ap.add_argument(
        "--show_fig_title", action="store_true",
        help="Writes SYSTEM as the figure's supertitle",
    )
    return ap

def get_raw_results(system : str) -> dict[str, float | str | NDArray]:
    """
    Loads the raw results for a system
    """
    system = system.upper()
    assert system in SYSTMS

    out = loadmat(f"./{system}_results", squeeze_me=True, simplify_cells=True)

    # Remove unnecessary entries
    for k in [k for k in out if (k.startswith("__") and k.endswith("__"))]: del out[k]
    for k in [k for k in out if not (
        k.endswith("ID") or k.endswith("EVALS") or k.endswith("REGS") or k.endswith("Oks")
        or k.startswith("Pu") or k.startswith("Qu")
        or k.startswith("Vm") or k.startswith("Va") or k.startswith("Ps")
    )]:
        del out[k]

    return out

def summarize_results(
    system : str = SYSTMS[0],
    slack_scale : int | None = None,
    marker_face_color : str = "blue",
    show_fig_title : bool = False,
) -> None:
    """
    Summarizes and plots the results for a specified system
    """
    system = system.upper()
    assert system in SYSTMS

    raw = get_raw_results(system)
    print(f"Summarizing results for {system}")

    dif_p = np.linalg.norm(raw['Pu'] - raw['Pu_snap'], ord=2, axis=1)
    dif_q = np.linalg.norm(raw['Qu'] - raw['Qu_snap'], ord=2, axis=1)

    FNAME = f"sum_{system}"
    fig = plt.figure(figsize=(6.03, 7.7), dpi=100)
    gsp = fig.add_gridspec(nrows=3, ncols=1)
    axs_p = fig.add_subplot(gsp[0, 0])
    axs_q = fig.add_subplot(gsp[1, 0])
    axs_s = fig.add_subplot(gsp[2, 0])

    # Active supply injections
    with plt.style.context("seaborn-pastel"):
        axs_p.plot(
            np.arange(raw['PUREGS'].shape[0]),
            dif_p,
            linestyle="None",
            marker="X", markersize=7.0,
            markeredgecolor="none",
            markerfacecolor=marker_face_color,
        )
    axs_p.yaxis.set_major_locator(plt.MaxNLocator(8))
    axs_p.grid(visible=True, axis="y", which="both")
    axs_p.set_ylim(bottom=None, top=None)
    axs_p.set_ylabel(r"$\alpha \left(\mu\right)$ [p.u.]", fontsize=8)
    axs_p.set_xticks(
        np.arange(raw['PUREGS'].shape[0]),
        labels=[f"({p:.1E},\n {q:.1E})" for p, q in zip(raw['PUREGS'], raw['QUREGS'])],
        fontsize=7,
    )
    for l in axs_p.get_yaxis().get_ticklabels(): l.set_fontsize(7)

    # Reactive supply injections
    with plt.style.context("seaborn-pastel"):
        axs_q.plot(
            np.arange(raw['QUREGS'].shape[0]),
            dif_q,
            linestyle="None",
            marker="d", markersize=7.0,
            markeredgecolor="none",
            markerfacecolor=marker_face_color,
        )
    axs_q.yaxis.set_major_locator(plt.MaxNLocator(8))
    axs_q.grid(visible=True, axis="y", which="both")
    axs_q.set_ylim(bottom=None, top=None)
    axs_q.set_ylabel(r"$\beta \left(\mu\right)$ [p.u.]", fontsize=8)
    axs_q.set_xticks(
        np.arange(raw['PUREGS'].shape[0]),
        labels=[f"({p:.1E},\n {q:.1E})" for p, q in zip(raw['PUREGS'], raw['QUREGS'])],
        fontsize=7,
    )
    for l in axs_q.get_yaxis().get_ticklabels(): l.set_fontsize(7)

    # Distributed slack
    yunit = r"$\times 10^{" + f"{-slack_scale}" + r"}$ " if (slack_scale and slack_scale > 0.) else ""
    with plt.style.context("seaborn-pastel"):
        axs_s.plot(
            np.arange(raw['PUREGS'].shape[0]),
            (10**(slack_scale or 0.0)) * raw['Ps'],
            linestyle="None",
            marker="^", markersize=7.0,
            markeredgecolor="none",
            markerfacecolor=marker_face_color,
        )
    axs_s.yaxis.set_major_locator(plt.MaxNLocator(8))
    axs_s.grid(visible=True, axis="y", which="both")
    axs_s.set_ylim(bottom=None, top=None)
    axs_s.set_xlabel(
        r"Supply regularization, $\left(\mu_{\mathrm{p}}, \mu_{\mathrm{q}}\right)$",
        fontsize=8,
    )
    axs_s.set_ylabel(
        r"Distibuted slack, $\kappa \left(\mu\right)$ [" + yunit + r"p.u.]",
        fontsize=8,
    )
    axs_s.set_xticks(
        np.arange(raw['PUREGS'].shape[0]),
        labels=[f"({p:.1E},\n {q:.1E})" for p, q in zip(raw['PUREGS'], raw['QUREGS'])],
        fontsize=7,
    )
    for l in axs_s.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    
    # Figure title
    if show_fig_title:
        fig.suptitle(system, y=0.963, fontsize=8, fontfamily="monospace")

    # Outro
    fig.align_ylabels(axs=[axs_p, axs_q, axs_s])
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.19)
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots saved to '{FNAME}.pdf'")

if __name__ == "__main__":
    args = get_argparser().parse_args()
    summarize_results(
        system=args.system,
        slack_scale=args.slack_scale,
        marker_face_color=args.marker_face_color,
        show_fig_title=args.show_fig_title,
    )

