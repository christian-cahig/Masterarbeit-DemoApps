"""
A script summarizing the results of probabilistic power flow analysis

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
from seaborn import histplot, kdeplot

__author__ = "Christian Cahig"
__version__ = "0.1.0"
__all__ = []

RawResults = dict[str, float | str | NDArray]
SYSTEMS = ["AEP57"]
XED_TYPES = ["uu", "ur", "ru", "rr"]

def get_argparser() -> ArgumentParser:
    """Loads the argument parser"""
    ap = ArgumentParser(description="A script for summarizing the results")
    ap.add_argument(
        "--system", type=str, required=True, choices=SYSTEMS,
        help="System whose results are to be summarized",
    )
    ap.add_argument(
        "--num_bins", type=int, required=False, default=None,
        help="Number of bins for the histograms",
    )
    return ap

def get_raw_results(
    system : str = SYSTEMS[0],
    drop_infeasible : bool = True,
) -> RawResults:
    """Loads the raw results for a system"""
    system = system.upper()
    out = loadmat(f"./{system}_results", squeeze_me=True, simplify_cells=True)

    # Remove unnecessary entries
    for k in [k for k in out if (k.startswith("__") and k.endswith("__"))]: del out[k]
    for k in [k for k in out if not (
        k.endswith("ID") or k.endswith("MVA") or k.endswith("BUSES")
        or k.endswith("Pds") or k.endswith("Qds")
        or k.startswith("uu") or k.startswith("ru")
        or k.startswith("ur") or k.startswith("rr")
    )]: del out[k]

    # Remove infeasible entries
    if drop_infeasible:
        for t in XED_TYPES:
            idxs = np.bitwise_and(out[f"{t}Feasi"].astype(bool), out[f"{t}APFOk"].astype(bool))
            for e in ["Feasi", "APFOk", "Pus", "Qus", "Vms", "Vas", "Ps"]:
                out[f"{t}{e}"] = out[f"{t}{e}"][idxs]

    return out

def summarize_demand_draws(
    raw : RawResults,
    num_bins : int | str = "auto",
) -> plt.Figure:
    """Plot the distributions of the total active and of the total demand draws"""
    FNAME = f"PdQd_{raw['SNAPSHOT_ID']}"
    fig = plt.figure(figsize=(5.99, 7.0), dpi=100)
    gsp = fig.add_gridspec(nrows=2, ncols=1)
    axs_p = fig.add_subplot(gsp[0, 0])
    axs_q = fig.add_subplot(gsp[1, 0])

    # Active demand draws
    x = raw['meanPds'].sum() * raw['multPds']
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x, kde=True, bins=num_bins, color="lightcoral", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs_p,
        )
    axs_p.axvline(x.mean(), color="lightcoral", linestyle="--", linewidth=1.7, label="Mean")
    axs_p.grid(visible=True, axis="y", which="major")
    axs_p.set_xlabel("Expected total active demand draw [p.u]", fontsize=8)
    axs_p.set_ylabel("Count", fontsize=8)
    for l in axs_p.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_p.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    axs_p.legend(loc="upper right", fontsize=8, shadow=True)

    # Reactive demand draws
    x = raw['meanQds'].sum() * raw['multQds']
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x, kde=True, bins=num_bins, color="sienna", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs_q,
        )
    axs_q.axvline(x.mean(), color="sienna", linestyle="--", linewidth=1.7, label="Mean")
    axs_q.grid(visible=True, axis="y", which="major")
    axs_q.set_xlabel("Expected total reactive demand draw [p.u.]", fontsize=8)
    axs_q.set_ylabel("Count", fontsize=8)
    for l in axs_q.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_q.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    axs_q.legend(loc="upper right", fontsize=8, shadow=True)

    # Outro
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.17)
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots saved to '{FNAME}.pdf'")
    return fig

def summarize_supply_injections(
    raw : RawResults,
    num_bins : int | str = "auto",
    xed_type : str = XED_TYPES[0],
) -> plt.Figure:
    """Plot the distributions of the total active and of the total supply injections"""
    xed_type = xed_type.lower()
    assert xed_type in XED_TYPES

    FNAME = f"PuQu-{xed_type}_{raw['SNAPSHOT_ID']}"
    fig = plt.figure(figsize=(5.99, 7.0), dpi=100)
    gsp = fig.add_gridspec(nrows=2, ncols=1)
    axs_p = fig.add_subplot(gsp[0, 0])
    axs_q = fig.add_subplot(gsp[1, 0])

    # Active supply injections
    x = raw[f"{xed_type}Pus"].sum(axis=-1)
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x, kde=True, bins=num_bins, color="dodgerblue", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs_p,
        )
    axs_p.axvline(x.mean(), color="dodgerblue", linestyle="--", linewidth=1.7, label="Mean")
    axs_p.grid(visible=True, axis="y", which="major")
    axs_p.set_xlabel("Total ancitipated active supply injection [p.u.]", fontsize=8)
    axs_p.set_ylabel("Count", fontsize=8)
    for l in axs_p.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_p.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    axs_p.legend(loc="upper right", fontsize=8, shadow=True)

    # Reactive supply injections
    x = raw[f"{xed_type}Qus"].sum(axis=-1)
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x, kde=True, bins=num_bins, color="forestgreen", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs_q,
        )
    axs_q.axvline(x.mean(), color="forestgreen", linestyle="--", linewidth=1.7, label="Mean")
    axs_q.grid(visible=True, axis="y", which="major")
    axs_q.set_xlabel("Total ancitipated reactive supply injection [p.u.]", fontsize=8)
    axs_q.set_ylabel("Count", fontsize=8)
    for l in axs_q.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_q.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    axs_q.legend(loc="upper right", fontsize=8, shadow=True)

    # Outro
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.17)
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots saved to '{FNAME}.pdf'")
    return fig

def summarize_distributed_slack(
    raw : RawResults,
    num_bins : int | str = "auto",
    xed_type : str = XED_TYPES[0],
) -> plt.Figure:
    """Plot the distribution of the distributed slack variable"""
    xed_type = xed_type.lower()
    assert xed_type in XED_TYPES

    FNAME = f"Ps-{xed_type}_{raw['SNAPSHOT_ID']}"
    fig = plt.figure(figsize=(5.99, 3.5), dpi=100)
    gsp = fig.add_gridspec(nrows=1, ncols=1)
    axs = fig.add_subplot(gsp[0, 0])

    # Distributed slack
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=raw[f"{xed_type}Ps"], kde=True, bins=num_bins, color="teal", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs,
        )
    axs.axvline(
        raw[f"{xed_type}Ps"].mean(),
        color="teal", linestyle="--", linewidth=1.7,
        label="Mean"
    )
    axs.grid(visible=True, axis="y", which="major")
    axs.set_xlabel("Anticipated distributed slack [p.u.]", fontsize=8)
    axs.set_ylabel("Count", fontsize=8)
    for l in axs.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    axs.legend(loc="upper right", fontsize=8, shadow=True)

    # Outro
    fig.tight_layout()
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots saved to '{FNAME}.pdf'")
    return fig

if __name__ == "__main__":
    args = get_argparser().parse_args()
    raw = get_raw_results(system=args.system, drop_infeasible=True)
    summarize_demand_draws(raw, num_bins=args.num_bins if args.num_bins else "auto")
    for xed_type in XED_TYPES: summarize_supply_injections(
        raw,
        num_bins=args.num_bins if args.num_bins else "auto",
        xed_type=xed_type
    )
    for xed_type in XED_TYPES: summarize_distributed_slack(
        raw,
        num_bins=args.num_bins if args.num_bins else "auto",
        xed_type=xed_type
    )
