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
        "system", type=str, choices=SYSTEMS,
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

def summarize_voltage_mags(
    raw : RawResults,
    num_bins : int | str = "auto",
    xed_type : str = XED_TYPES[0],
) -> plt.Figure:
    """Plot the distribution of the bus-voltage magnitudes"""
    xed_type = xed_type.lower()
    assert xed_type in XED_TYPES

    x0 = raw[f"{xed_type}Vms"].mean(axis=0)
    idx1, idx2 = x0.argmin(), x0.argmax()

    FNAME = f"Vm-{xed_type}_{raw['SNAPSHOT_ID']}"
    fig = plt.figure(figsize=(5.99, 10), dpi=100)
    gsp = fig.add_gridspec(nrows=3, ncols=1)
    axs0 = fig.add_subplot(gsp[0, 0])
    axs1 = fig.add_subplot(gsp[1, 0])
    axs2 = fig.add_subplot(gsp[2, 0])

    # Over all buses
    x0 = raw[f"{xed_type}Vms"].flatten()
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x0, kde=True, bins=num_bins, color="seagreen", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs0,
        )
    axs0.set_xlabel("Anticipated bus-voltage magnitudes [p.u.]", fontsize=8)

    # Over the bus with the least average voltage magnitude
    x1 = raw[f"{xed_type}Vms"][:, idx1]
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x1, kde=True, bins=num_bins, color="seagreen", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs1,
        )
    axs1.set_xlabel(f"Anticipated voltage magnitude [p.u.] at bus {idx1+1}", fontsize=8)

    # Over the bus with the largest average voltage magnitude
    x2 = raw[f"{xed_type}Vms"][:, idx2]
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x2, kde=True, bins=num_bins, color="seagreen", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs2,
        )
    axs2.set_xlabel(f"Anticipated voltage magnitude [p.u.] at bus {idx2+1}", fontsize=8)

    # Visuals
    for axs, x, in zip([axs0, axs1, axs2], [x0, x1, x2]):
        axs.axvline(x.mean(), color="seagreen", linestyle="--", linewidth=1.7, label="Mean")
        axs.grid(visible=True, axis="y", which="major")
        axs.set_xlim(left=0.87, right=1.13)
        axs.set_ylabel("Count", fontsize=8)
        for l in axs.get_xaxis().get_ticklabels(): l.set_fontsize(7)
        for l in axs.get_yaxis().get_ticklabels(): l.set_fontsize(7)
        axs.legend(loc="upper right", fontsize=8, shadow=True)

    # Outro
    fig.align_ylabels(axs=[axs0, axs1, axs2])
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.17)
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots saved to '{FNAME}.pdf'")
    return fig

def summarize_voltage_angs(
    raw : RawResults,
    num_bins : int | str = "auto",
    xed_type : str = XED_TYPES[0],
) -> plt.Figure:
    """Plot the distribution of the branch-phase-angle differences"""
    xed_type = xed_type.lower()
    assert xed_type in XED_TYPES

    x0 = raw[f"{xed_type}Vas"].mean(axis=0)
    idx1, idx2 = x0.argmin(), x0.argmax()

    FNAME = f"Va-{xed_type}_{raw['SNAPSHOT_ID']}"
    fig = plt.figure(figsize=(5.99, 10), dpi=100)
    gsp = fig.add_gridspec(nrows=3, ncols=1)
    axs0 = fig.add_subplot(gsp[0, 0])
    axs1 = fig.add_subplot(gsp[1, 0])
    axs2 = fig.add_subplot(gsp[2, 0])

    # Over all branches
    x0 = raw[f"{xed_type}Vas"].flatten()
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x0, kde=True, bins=num_bins, color="cornflowerblue", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs0,
        )
    axs0.set_xlim(left=-0.067, right=0.137)
    axs0.set_xlabel("Anticipated branch-phase-angle differences [°]", fontsize=8)

    # Over the branch with the least average voltage-phase-angle difference
    x1 = raw[f"{xed_type}Vas"][:, idx1]
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x1, kde=True, bins=num_bins, color="cornflowerblue", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs1,
        )
    axs1.set_xlabel(
        ("Anticipated phase-angle difference [°] at branch connecting "
        f"buses {raw['FBUSES'][idx1]} and {raw['TBUSES'][idx1]}"),
        fontsize=8
    )

    # Over the branch with the largest average voltage-phase-angle difference
    x2 = raw[f"{xed_type}Vas"][:, idx2]
    with plt.style.context("seaborn-pastel"):
        histplot(
            x=x2, kde=True, bins=num_bins, color="cornflowerblue", fill=False,
            line_kws={"linestyle" : "-.", "linewidth" : 1.0, "label" : "KDE"},
            ax=axs2,
        )
    axs2.set_xlabel(
        ("Anticipated phase-angle difference [°] at branch connecting "
        f"buses {raw['FBUSES'][idx2]} and {raw['TBUSES'][idx2]}"),
        fontsize=8
    )

    # Visuals
    for axs, x, in zip([axs0, axs1, axs2], [x0, x1, x2]):
        axs.axvline(x.mean(), color="cornflowerblue", linestyle="--", linewidth=1.7, label="Mean")
        axs.grid(visible=True, axis="y", which="major")
        axs.set_xlim(left=None, right=None)
        axs.set_ylabel("Count", fontsize=8)
        for l in axs.get_xaxis().get_ticklabels(): l.set_fontsize(7)
        for l in axs.get_yaxis().get_ticklabels(): l.set_fontsize(7)
        axs.legend(loc="upper right", fontsize=8, shadow=True)

    # Outro
    fig.align_ylabels(axs=[axs0, axs1, axs2])
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.17)
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots saved to '{FNAME}.pdf'")
    return fig

if __name__ == "__main__":
    args = get_argparser().parse_args()
    n = args.num_bins if args.num_bins else "auto"

    # Results data
    raw = get_raw_results(system=args.system, drop_infeasible=True)

    # Summarize demand draws
    print("Averages of total active demand")
    print(f"{' '*3}{(raw['meanPds'].sum() * raw['multPds']).mean()}")
    print("Standard deviations of total active demand")
    print(f"{' '*3}{(raw['meanPds'].sum() * raw['multPds']).std()}")
    print("Averages of total reactive demand")
    print(f"{' '*3}{(raw['meanQds'].sum() * raw['multQds']).mean()}")
    print("Standard deviations of total reactive demand")
    print(f"{' '*3}{(raw['meanQds'].sum() * raw['multQds']).std()}")
    summarize_demand_draws(raw, num_bins=n)

    # Summarize anticipated supply injections
    print("Averages of total active supply")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Pus"].sum(axis=-1).mean(),
        sep=": "
    )
    print("Standard deviations of total active supply")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Pus"].sum(axis=-1).std(),
        sep=": "
    )
    print("Averages of total reactive supply")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Qus"].sum(axis=-1).mean(),
        sep=": "
    )
    print("Standard deviations of total reactive supply")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Qus"].sum(axis=-1).std(),
        sep=": "
    )
    for t in XED_TYPES: summarize_supply_injections(raw, num_bins=n, xed_type=t)

    # Summarize anticipated bus voltage magnitudes
    print("Averages over all bus-voltage magnitudes")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Vms"].flatten().mean(),
        sep=": "
    )
    print("Standard deviations over all bus-voltage magnitudes")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Vms"].flatten().std(),
        sep=": "
    )
    print("Averages over bus with least average voltage magnitude")
    for xed_type in XED_TYPES:
        x = raw[f"{xed_type}Vms"].mean(axis=0)
        idx = x.argmin()
        print(
            f"{' '*3}{xed_type} @ {idx+1}",
            raw[f"{xed_type}Vms"][:, idx].mean(),
            sep=": "
        )
    print(f"Standard deviations over bus with least average voltage magnitude")
    for xed_type in XED_TYPES:
        x = raw[f"{xed_type}Vms"].mean(axis=0)
        idx = x.argmin()
        print(
            f"{' '*3}{xed_type} @ {idx+1}",
            raw[f"{xed_type}Vms"][:, idx].std(),
            sep=": "
        )
    print("Averages over bus with largest average voltage magnitude")
    for xed_type in XED_TYPES:
        x = raw[f"{xed_type}Vms"].mean(axis=0)
        idx = x.argmax()
        print(
            f"{' '*3}{xed_type} @ {idx+1}",
            raw[f"{xed_type}Vms"][:, idx].mean(),
            sep=": "
        )
    print(f"Standard deviations over bus with largest average voltage magnitude")
    for xed_type in XED_TYPES:
        x = raw[f"{xed_type}Vms"].mean(axis=0)
        idx = x.argmax()
        print(
            f"{' '*3}{xed_type} @ {idx+1}",
            raw[f"{xed_type}Vms"][:, idx].std(),
            sep=": "
        )
    for t in XED_TYPES: summarize_voltage_mags(raw, num_bins=n, xed_type=t)

    # Summarize anticipated branch phase-angle differences
    print("Averages over all branch-phase-angle differences")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Vas"].flatten().mean(),
        sep=": "
    )
    print("Standard deviations over all branch-phase-angle differences")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Vas"].flatten().std(),
        sep=": "
    )
    print("Averages over branch with least average phase-angle difference")
    for xed_type in XED_TYPES:
        x = raw[f"{xed_type}Vas"].mean(axis=0)
        idx = x.argmin()
        print(
            f"{' '*3}{xed_type} @ ({raw['FBUSES'][idx]} -> {raw['TBUSES'][idx]})",
            raw[f"{xed_type}Vas"][:, idx].mean(),
            sep=": "
        )
    print(f"Standard deviations over branch with least average phase-angle difference")
    for xed_type in XED_TYPES:
        x = raw[f"{xed_type}Vas"].mean(axis=0)
        idx = x.argmin()
        print(
            f"{' '*3}{xed_type} @ ({raw['FBUSES'][idx]} -> {raw['TBUSES'][idx]})",
            raw[f"{xed_type}Vas"][:, idx].std(),
            sep=": "
        )
    print("Averages over branch with largest average phase-angle difference")
    for xed_type in XED_TYPES:
        x = raw[f"{xed_type}Vas"].mean(axis=0)
        idx = x.argmax()
        print(
            f"{' '*3}{xed_type} @ ({raw['FBUSES'][idx]} -> {raw['TBUSES'][idx]})",
            raw[f"{xed_type}Vas"][:, idx].mean(),
            sep=": "
        )
    print(f"Standard deviations over branch with largest average phase-angle difference")
    for xed_type in XED_TYPES:
        x = raw[f"{xed_type}Vas"].mean(axis=0)
        idx = x.argmax()
        print(
            f"{' '*3}{xed_type} @ ({raw['FBUSES'][idx]} -> {raw['TBUSES'][idx]})",
            raw[f"{xed_type}Vas"][:, idx].std(),
            sep=": "
        )
    for t in XED_TYPES: summarize_voltage_angs(raw, num_bins=n, xed_type=t)

    # Summarize anticipated distributed slacks
    print("Averages of anticipated distributed slack")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Ps"].mean(),
        sep=": "
    )
    print("Standard deviations of anticipated distributed slack")
    for xed_type in XED_TYPES: print(
        f"{' '*3}{xed_type}",
        raw[f"{xed_type}Ps"].std(),
        sep=": "
    )
    for t in XED_TYPES: summarize_distributed_slack(raw, num_bins=n, xed_type=t)
