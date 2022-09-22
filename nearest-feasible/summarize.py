"""
A script summarizing the results of using APF equations to find nearest power flow feasible points

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
        "--num_instances", type=int, required=False, default=None,
        help="Number of instances to show",
    )
    ap.add_argument(
        "--slk_scale", type=int, required=False, default=None,
        help="Scales the distributed slacks by 10^(SLK_SCALE)",
    )
    ap.add_argument(
        "--mag_scale", type=int, required=False, default=None,
        help="Scales the differences in voltage magnitude profiles by 10^(MAG_SCALE)",
    )
    ap.add_argument(
        "--ang_scale", type=int, required=False, default=None,
        help="Scales the differences in voltage-angle-difference profiles by 10^(ANG_SCALE)",
    )
    return ap

def get_raw_result(
    system : str = SYSTMS[0],
    drop_infeasible : bool = False,
    num_instances : int | None = None,
    slk_scale : int | None = None,
    mag_scale : int | None = None,
    ang_scale : int | None = None,
) -> RawResultDict:
    """
    Loads the raw results
    """
    system = system.upper()
    assert system in SYSTMS

    out = loadmat(f"./{system}_results.mat", squeeze_me=True, simplify_cells=True)

    # Remove unnecessary entries
    for k in [k for k in out if (k.startswith("__") and k.endswith("__"))]: del out[k]

    # Convert to booleans
    out['apfOks'] = out['apfOks'].astype(bool)

    # Drop power-flow-infeasible instances
    if drop_infeasible:
        for k in [k for k in out if (
            k.endswith("Pus") or k.endswith("Qus")
            or k.endswith("Vms") or k.endswith("Vas")
            or k.endswith("Vbs") or k.endswith("Sks")
        )]:
            out[k] = out[k][out['apfOks']]
        out['apfOks'] = out['apfOks'][out['apfOks']]

    # Scale distributed slacks
    slk_scale = slk_scale or 0.
    out['apfSks'] = (10**(-slk_scale)) * out['apfSks']

    # Scale differences in voltage magnitude profiles
    mag_scale = mag_scale or 0.
    out['difVms'] = (10**(-mag_scale)) * out['difVms']

    # Scale differences in voltage-angle-difference profiles
    ang_scale = ang_scale or 0.
    out['difVas'] = (10**(-ang_scale)) * out['difVas']

    # Keep only the specified number of instances
    if num_instances:
        num_instances = max(10, min(num_instances, out['apfOks'].shape[0]))
        idxs = np.random.choice(
            np.arange(out['apfOks'].shape[0]),
            size=(num_instances,),
            replace=False,
        )
        for k in [k for k in out if (
            k.endswith("Oks")
            or k.endswith("Pus") or k.endswith("Qus")
            or k.endswith("Vms") or k.endswith("Vas")
            or k.endswith("Vbs") or k.endswith("Vds") or k.endswith("Sks")
        )]:
            out[k] = out[k][idxs]

    # Add system name
    out['SYS'] = system

    return out

def plot_summary(args : Namespace) -> plt.Figure:
    """
    Plots and saves the results
    """
    raw = get_raw_result(
        system=args.system, drop_infeasible=True,
        slk_scale=args.slk_scale,
        mag_scale=args.mag_scale,
        ang_scale=args.ang_scale,
        num_instances=args.num_instances,
    )
    print(f"Summarizing results for {raw['SYS']}")

    FNAME = f"approx-sum_{raw['SYS']}"
    fig = plt.figure(figsize=(6.015, 7.25), dpi=100)
    gsp = fig.add_gridspec(nrows=3, ncols=1)
    axs_s = fig.add_subplot(gsp[0, 0])
    axs_v = fig.add_subplot(gsp[1, 0])
    axs_a = fig.add_subplot(gsp[2, 0])

    # Distributed slacks
    _ave = raw['apfSks'].mean()
    _std = raw['apfSks'].std()
    yunit = r"$\times 10^{" + f"{args.slk_scale}" + r"}$ " if (
        args.slk_scale and (args.slk_scale != 0.)
    ) else ""
    with plt.style.context("seaborn-pastel"):
        axs_s.axhline(
            y=_ave,
            linestyle="--", linewidth=1.0,
            color="darkolivegreen",
        )
        axs_s.plot(
            np.arange(1, raw['apfOks'].shape[0] + 1),
            raw['apfSks'],
            linestyle="None",
            marker="1", markersize=8.0, markeredgewidth=0.7,
            markerfacecolor="none",
            markeredgecolor="darkolivegreen",
        )
    axs_s.grid(visible=True, axis="y", which="both")
    axs_s.set_xlim(left=None, right=None)
    axs_s.set_ylim(bottom=None, top=None)
    axs_s.set_ylabel(r"Distributed slack, $\kappa$ [" + yunit + r"p.u.]", fontsize=8)
    axs_s.set_xticklabels([])
    axs_s.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    for l in axs_s.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    print(f"{' ' * 3}Distributed slack [p.u.]: ", end="")
    print(f"{_ave} ± {_std:E}")

    # Voltage deviation indices
    _ave = raw['difVms'].mean()
    _std = raw['difVms'].std()
    yunit = r"$\times 10^{" + f"{args.mag_scale}" + r"}$ " if (
        args.mag_scale and (args.mag_scale != 0.)
    ) else ""
    with plt.style.context("seaborn-pastel"):
        axs_v.axhline(
            y=_ave,
            linestyle="--", linewidth=1.0,
            color="sandybrown",
        )
        axs_v.plot(
            np.arange(1, raw['apfOks'].shape[0] + 1),
            raw['difVms'],
            linestyle="None",
            marker="*", markersize=7.0, markeredgewidth=0.7,
            markerfacecolor="none",
            markeredgecolor="sandybrown",
        )
    axs_v.grid(visible=True, axis="y", which="both")
    axs_v.set_xlim(left=None, right=None)
    axs_v.set_ylim(bottom=None, top=None)
    axs_v.set_ylabel(r"$\epsilon_{\mathrm{v}}$ [" + yunit + r"p.u.]", fontsize=8)
    axs_v.set_xticklabels([])
    axs_v.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    for l in axs_v.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    print(f"{' ' * 3}Bus voltage magnitudes [p.u.]: ", end="")
    print(f"{_ave} ± {_std:E}")

    # Differences in voltage-angle-difference profiles
    _ave = raw['difVas'].mean()
    _std = raw['difVas'].std()
    yunit = r"$\times 10^{" + f"{args.ang_scale}" + r"}$ " if (
        args.ang_scale and (args.ang_scale != 0.)
    ) else ""
    with plt.style.context("seaborn-pastel"):
        axs_a.axhline(
            y=_ave,
            linestyle="--", linewidth=1.0,
            color="mediumpurple",
        )
        axs_a.plot(
            np.arange(1, raw['apfOks'].shape[0] + 1),
            raw['difVas'],
            linestyle="None",
            marker="o", markersize=5.0, markeredgewidth=0.7,
            markerfacecolor="none",
            markeredgecolor="mediumpurple",
        )
    axs_a.grid(visible=True, axis="y", which="both")
    axs_a.set_xlim(left=None, right=None)
    axs_a.set_ylim(bottom=None, top=None)
    axs_a.set_xlabel("Augmented DCOPF instance", fontsize=8)
    axs_a.set_ylabel(r"$\epsilon_{\mathrm{a}}$ [" + yunit + r"rad.]", fontsize=8)
    axs_a.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    for l in axs_a.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_a.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    print(f"{' ' * 3}Branch voltage phase angle differences [degrees]: ", end="")
    print(f"{np.rad2deg(_ave)} ± {np.rad2deg(_std):E}")

    # Outro
    fig.align_ylabels(axs=[axs_s, axs_v, axs_a])
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.05)
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots saved to '{FNAME}.pdf'")

    return fig

if __name__ == "__main__":
    plot_summary(get_argparser().parse_args())
