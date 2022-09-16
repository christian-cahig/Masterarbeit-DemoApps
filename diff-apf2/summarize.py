"""
A script summarizing the results of differentiating through the APF equations

Masterarbeit-DemApp
Copyright (C) 2022 - present, Christian Cahig
https://github.com/christian-cahig/Masterarbeit-DemApp

This file is part of Masterarbeit-DemApp, and is covered by the CC-BY-4.0 License.
See `LICENSE` file in the source directory for details.
"""

from argparse import ArgumentParser, Namespace
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

import numpy as np
from numpy.typing import NDArray
from scipy.io import loadmat
from scipy.sparse import spmatrix as SpMatrix
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

__author__ = "Christian Cahig"
__version__ = "0.1.0"
__all__ = []

SYSTMS = ["ACT500", "RTE1888", "POL3375"]
RawResultDict = dict[str, str | int | float | NDArray | SpMatrix]

def get_argparser() -> ArgumentParser:
    """
    Returns the argument parser
    """
    ap = ArgumentParser(description="A script for summarizing the results")
    ap.add_argument(
        "system", type=str, choices=SYSTMS,
        help="System whose results are to be summarized",
    )
    ap.add_argument(
        "--type", type=str, required=False, default="dPs",
        choices={"dPs", "dF", "tdPs", "tdF"},
        help="Result type to be summarized",
    )
    ap.add_argument(
        "--pscale", type=int, required=False, default=None,
        help="Scales quantities pertaining to active power by 10^(PSCALE)",
    )
    ap.add_argument(
        "--qscale", type=int, required=False, default=None,
        help="Scales quantities pertaining to reactive power by 10^(QSCALE)",
    )
    ap.add_argument(
        "--pcolor", type=str, required=False, default="coral",
        help="Color for plots pertaining to active power",
    )
    ap.add_argument(
        "--qcolor", type=str, required=False, default="mediumpurple",
        help="Color for plots pertaining to reactive power",
    )
    ap.add_argument(
        "--pticks", type=int, required=False, default=5,
        help="Number of y-axis ticks in plots pertaining to active power",
    )
    ap.add_argument(
        "--qticks", type=int, required=False, default=5,
        help="Number of y-axis ticks in plots pertaining to reactive power",
    )
    return ap

def get_raw_result(system : str = SYSTMS[0]) -> RawResultDict:
    """
    Loads the raw results for a specified system
    """
    system = system.upper()
    assert system in SYSTMS

    out = loadmat(f"./{system}_presults.mat", squeeze_me=True, simplify_cells=True)

    # Remove unnecessary entries
    for k in [k for k in out if (k.startswith("__") and k.endswith("__"))]: del out[k]

    # Fix vectors in sparse-matrix format
    out['dPsdC'] = out['dPsdC'].data

    # Add system name
    out['SYS'] = system

    return out

def plot_dPs(
    raw : RawResultDict,
    pscale : int | None = None,
    qscale : int | None = None,
    pcolor : str = "coral",
    qcolor : str = "mediumpurple",
    pticks : int | None = None,
    qticks : int | None = None,
) -> plt.Figure:
    """
    Plots and saves results concerning the gradient of the distributed slack
    """
    pscale, qscale = pscale or 0., qscale or 0.

    FNAME = f"diff-dPs_{raw['SYS']}"
    fig = plt.figure(figsize=(6.01, 4.70), dpi=100)
    gsp = fig.add_gridspec(nrows=2, ncols=1)
    axs_p = fig.add_subplot(gsp[0, 0])
    axs_q = fig.add_subplot(gsp[1, 0])
    x = np.arange(1, raw['N_u'] + 1)

    # Mismatches of gradients w.r.t. anticipated active supply injections
    y = 10**(-pscale) * (raw['dPsdC'][:raw['N_u']] - raw['dPsdC_'][:raw['N_u']])
    print(f"{' ' * 3}Min. dPsdPu (scaled): {y.min()}")
    print(f"{' ' * 3}Max. dPsdPu (scaled): {y.max()}")
    print(f"{' ' * 3}Ave. dPsdPu (scaled): {y.mean()}")
    yunit = r"$\left[ \times 10^{" + f"{pscale}" + r"} \right]$" if (
        pscale and (pscale != 0.)
    ) else ""
    with plt.style.context("seaborn-pastel"):
        axs_p.plot(
            x, y,
            linestyle="-", linewidth=0.5,
            color=pcolor,
            marker="4", markersize=7.0, markeredgewidth=0.3,
            markerfacecolor="none",
            markeredgecolor=pcolor,
        )
    if pticks: axs_p.yaxis.set_major_locator(plt.MaxNLocator(pticks))
    axs_p.grid(visible=True, axis="y", which="both")
    axs_p.set_xlim(left=None, right=None)
    axs_p.set_ylabel(
        r"$\Delta\kappa_{\mathrm{p},\!i}$ " + yunit,
        fontsize=8,
    )
    axs_p.set_xticklabels([])
    for l in axs_p.get_yaxis().get_ticklabels(): l.set_fontsize(7)

    # Mismatches of gradient w.r.t. anticipated reactive supply injections
    y = 10**(-pscale) * (raw['dPsdC'][raw['N_u']:] - raw['dPsdC_'][raw['N_u']:])
    print(f"{' ' * 3}Min. dPsdQu (scaled): {y.min()}")
    print(f"{' ' * 3}Max. dPsdQu (scaled): {y.max()}")
    print(f"{' ' * 3}Ave. dPsdQu (scaled): {y.mean()}")
    yunit = r"$\left[ \times 10^{" + f"{qscale}" + r"} \right]$" if (
        qscale and (qscale != 0.)
    ) else ""
    with plt.style.context("seaborn-pastel"):
        axs_q.plot(
            x, y,
            linestyle="-", linewidth=0.5,
            color=qcolor,
            marker="o", markersize=4.0, markeredgewidth=0.3,
            markerfacecolor="none",
            markeredgecolor=qcolor,
        )
    if qticks: axs_q.yaxis.set_major_locator(plt.MaxNLocator(qticks))
    axs_q.grid(visible=True, axis="y", which="both")
    axs_q.set_xlim(left=None, right=None)
    axs_q.set_xlabel(r"Supply unit, $i$", fontsize=8)
    axs_q.set_ylabel(
        r"$\Delta\kappa_{\mathrm{q},\!i}$ " + yunit,
        fontsize=8,
    )
    for l in axs_q.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_q.get_yaxis().get_ticklabels(): l.set_fontsize(7)

    # Outro
    fig.align_ylabels(axs=[axs_p, axs_q])
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.05)
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}dPsdC plots saved to '{FNAME}.pdf'")
    return fig

def plot_dF(
    raw : RawResultDict,
    pscale : int | None = None,
    qscale : int | None = None,
    pcolor : str = "coral",
    qcolor : str = "mediumpurple",
    pticks : int | None = None,
    qticks : int | None = None,
) -> plt.Figure:
    """
    Plots and saves results concerning the gradient of the distributed slack
    """
    pscale, qscale = pscale or 0., qscale or 0.

    FNAME = f"diff-dF_{raw['SYS']}"
    fig = plt.figure(figsize=(6.01, 4.70), dpi=100)
    gsp = fig.add_gridspec(nrows=2, ncols=1)
    axs_p = fig.add_subplot(gsp[0, 0])
    axs_q = fig.add_subplot(gsp[1, 0])
    x = np.arange(1, raw['N_u'] + 1)

    # Mismatches of gradients w.r.t. anticipated active supply injections
    y = 10**(-pscale) * (raw['dFdC'][:raw['N_u']] - raw['dFdC_'][:raw['N_u']])
    print(f"{' ' * 3}Min. dFdPu (scaled): {y.min()}")
    print(f"{' ' * 3}Max. dFdPu (scaled): {y.max()}")
    print(f"{' ' * 3}Ave. dFdPu (scaled): {y.mean()}")
    yunit = r"$\left[ \times 10^{" + f"{pscale}" + r"} \right]$" if (
        pscale and (pscale != 0.)
    ) else ""
    with plt.style.context("seaborn-pastel"):
        axs_p.plot(
            x, y,
            linestyle="-", linewidth=0.5,
            color=pcolor,
            marker="4", markersize=7.0, markeredgewidth=0.3,
            markerfacecolor="none",
            markeredgecolor=pcolor,
        )
    if pticks: axs_p.yaxis.set_major_locator(plt.MaxNLocator(pticks))
    axs_p.grid(visible=True, axis="y", which="both")
    axs_p.set_xlim(left=None, right=None)
    axs_p.set_ylabel(
        r"$\Delta\ell_{\mathrm{p},\!i}$ " + yunit,
        fontsize=8,
    )
    axs_p.set_xticklabels([])
    for l in axs_p.get_yaxis().get_ticklabels(): l.set_fontsize(7)

    # Mismatches of gradient w.r.t. anticipated reactive supply injections
    y = 10**(-pscale) * (raw['dFdC'][raw['N_u']:] - raw['dFdC_'][raw['N_u']:])
    print(f"{' ' * 3}Min. dFdQu (scaled): {y.min()}")
    print(f"{' ' * 3}Max. dFdQu (scaled): {y.max()}")
    print(f"{' ' * 3}Ave. dFdQu (scaled): {y.mean()}")
    yunit = r"$\left[ \times 10^{" + f"{qscale}" + r"} \right]$" if (
        qscale and (qscale != 0.)
    ) else ""
    with plt.style.context("seaborn-pastel"):
        axs_q.plot(
            x, y,
            linestyle="-", linewidth=0.5,
            color=qcolor,
            marker="o", markersize=4.0, markeredgewidth=0.3,
            markerfacecolor="none",
            markeredgecolor=qcolor,
        )
    if qticks: axs_q.yaxis.set_major_locator(plt.MaxNLocator(qticks))
    axs_q.grid(visible=True, axis="y", which="both")
    axs_q.set_xlim(left=None, right=None)
    axs_q.set_xlabel(r"Supply unit, $i$", fontsize=8)
    axs_q.set_ylabel(
        r"$\Delta\ell_{\mathrm{q},\!i}$ " + yunit,
        fontsize=8,
    )
    for l in axs_q.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_q.get_yaxis().get_ticklabels(): l.set_fontsize(7)

    # Outro
    fig.align_ylabels(axs=[axs_p, axs_q])
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.05)
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}dFdC plots saved to '{FNAME}.pdf'")
    return fig

def report_run_time(
    raw : RawResultDict,
    type : str = "Ps",
) -> None:
    """
    Prints a summary of backpropagation run times
    """
    assert type in {"Ps", "F"}

    if type == "Ps":
        print(f"{' '*3}Mismatch: {np.linalg.norm(raw['dPsdC'] - raw['dPsdC_'], ord=2)}")
        print(f"{' '*3}CFD: {raw['dPsdCTimes_'].mean()} ± {raw['dPsdCTimes_'].std()} s")
        print(f"{' '*3}IFT: {raw['dPsdCTimeAve']} ± {raw['dPsdCTimeStd']} s")
    else:
        print(f"{' '*3}Mismatch: {np.linalg.norm(raw['DFDC'] - raw['dFdC'], ord=np.Inf)}")
        print(f"{' '*3}BCR: {raw['dFdCTimeAve']} ± {raw['dFdCTimeStd']} s")
        print(f"{' '*3}VJP: {raw['DFDCTimeAve']} ± {raw['DFDCTimeStd']} s")

if __name__ == "__main__":
    args = get_argparser().parse_args()
    raw = get_raw_result(args.system)
    print(f"Summarizing {raw['SYS']}")

    if args.type == "dPs":
        plot_dPs(raw, pscale=args.pscale, qscale=args.pscale,
                 pcolor=args.pcolor, qcolor=args.qcolor,
                 pticks=args.pticks, qticks=args.qticks)
    elif args.type == "dF":
        plot_dF(raw, pscale=args.pscale, qscale=args.pscale,
                 pcolor=args.pcolor, qcolor=args.qcolor,
                 pticks=args.pticks, qticks=args.qticks)
    elif args.type == "tdPs":
        report_run_time(raw, type="Ps")
    elif args.type == "tdF":
        report_run_time(raw, type="F")
    else:
        raise NotImplementedError
