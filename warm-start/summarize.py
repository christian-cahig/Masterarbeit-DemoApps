"""
A script summarizing the results of using APF as a warm-starter for OPF routines

Masterarbeit-DemApp
Copyright (C) 2022 - present, Christian Cahig
https://github.com/christian-cahig/Masterarbeit-DemApp

This file is part of Masterarbeit-DemApp, and is covered by the CC-BY-4.0 License.
See `LICENSE` file in the source directory for details.
"""

from argparse import ArgumentParser
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

import numpy as np, pandas as pd
from numpy.typing import NDArray
from scipy.io import loadmat
import matplotlib.pyplot as plt

__author__ = "Christian Cahig"
__version__ = "0.1.0"
__all__ = []

SYSTMS = ["ACT500", "RTE1888", "POL3375"]
RawResultDict = dict[str, str | int | float | NDArray]

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
        "--sigfig", type=int, required=False, default=6,
        help="Number of significant digits to use",
    )
    return ap

def _float2sci(num : float, digits : int = 3, format : str = "g") -> str:
    """
    Returns a string representation of the scientific notation of a given number

    This is essentially a copy of T.C. Proctor's StackOverflow answer:
    https://stackoverflow.com/a/59744605
    """
    e_float = "{0:.{1:d}{}}".format(num, digits, format)
    if "e" not in e_float: return f"${e_float}$"
    mantissa, exponent = e_float.split("e")
    cleaned_exponent = exponent.strip("+")
    return f"${mantissa} \\times 10^{{{cleaned_exponent}}}$"

def get_raw_result(system : str = SYSTMS[0], drop_infeasible : bool = False) -> RawResultDict:
    """
    Loads the raw results concerning one system
    """
    system = system.upper()
    assert system in SYSTMS

    out = loadmat(f"./{system}_results.mat", squeeze_me=True, simplify_cells=True)

    # Remove unnecessary entries
    for k in [k for k in out if (k.startswith("__") and k.endswith("__"))]: del out[k]

    # Convert to booleans
    out['areFeasible'] = out['areFeasible'].astype(bool)

    # Convert to integers
    for k in [k for k in out if k.endswith("Iters")]: out[k] = out[k].astype(np.uint8)

    # Drop power-flow-infeasible instances
    if drop_infeasible:
        for k in [k for k in out if (
            k.endswith("Objs") or k.endswith("Cs")
            or k.endswith("Vms") or k.endswith("Vas")
        )]:
            out[k] = out[k][out['areFeasible']]
        out['areFeasible'] = out['areFeasible'][out['areFeasible']]

    # Add system name
    out['SYSTEM'] = system

    return out

def _plot_optim_injections(raw : RawResultDict, num_instances : int = 100) -> None:
    """
    Plots the differences in optimum supply injections
    """
    num_instances = min(raw['areFeasible'].sum(), max(4, num_instances))
    prefixs = ["uu", "ur", "ru", "rr"]
    wsmodes = np.arange(len(prefixs)) + 1
    supregs = [(0, 0), (0, 1), (1, 0), (1, 1)]

    FNAME = f"sum_{raw['SYSTEM']}"
    fig = plt.figure(figsize=(6.01, 7.5), dpi=100)
    gsp = fig.add_gridspec(nrows=1, ncols=1)
    axs = fig.add_subplot(gsp[0, 0])

    idxs_max = {raw[k].argmax() for k in raw if k.endswith("Cs")}
    idxs = np.array(list(set(range(raw['areFeasible'].sum())).difference(idxs_max)))
    idxs = np.concatenate((
        np.random.choice(idxs, size=(num_instances - len(idxs_max)), replace=False),
        np.array(list(idxs_max)),
    ), axis=0)
    del idxs_max

    # Scatter plot
    with plt.style.context("seaborn-whitegrid"):
        scatplot = axs.scatter(
            x=np.concatenate([raw[f"{p}diffCs"][idxs] for p in prefixs], axis=0),
            y=np.repeat(
                np.arange(1, idxs.shape[0] + 1)[:, None], len(prefixs), axis=1,
            ).T.flatten(),
            c=np.repeat(wsmodes, idxs.shape[0], axis=0),
            cmap="Spectral",
        )

    # Cosmetics
    axs.grid(visible=True, axis="x", which="major")
    axs.set_xlim(left=None, right=None)
    axs.set_ylim(bottom=None, top=None)
    axs.set_xlabel(
        r"Distance from snapshot-started optimum supply injections [p.u.], $\epsilon$",
        fontsize=8,
    )
    axs.set_ylabel(r"OPF instance index, $i$", fontsize=8)
    for l in axs.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs.get_yaxis().get_ticklabels(): l.set_fontsize(7)

    # Legend

    # Outro
    fig.tight_layout()
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots saved to '{FNAME}.pdf'")

def tabulate_result(raw : RawResultDict) -> pd.DataFrame:
    """
    Summarizes the results as a DataFrame
    """
    prefxs = ["uu", "ur", "ru", "rr"]
    regs = [(0, 0), (0, 1), (1, 0), (1, 1)]
    return pd.DataFrame(data={
        "preg" : np.array([r[0] for r in regs]),
        "qreg" : np.array([r[1] for r in regs]),
        "minC" : np.array([raw[f"{p}diffCs"].min() for p in prefxs]),
        "maxC" : np.array([raw[f"{p}diffCs"].max() for p in prefxs]),
        "minVm" : np.array([raw[f"{p}diffVms"].min() for p in prefxs]),
        "maxVm" : np.array([raw[f"{p}diffVms"].max() for p in prefxs]),
        "minVa" : np.array([raw[f"{p}diffVas"].min() for p in prefxs]),
        "maxVa" : np.array([raw[f"{p}diffVas"].max() for p in prefxs]),
        "minF" : np.array([(raw[f"{p}Objs"] - raw['snapObjs']).min() for p in prefxs]),
        "maxF" : np.array([(raw[f"{p}Objs"] - raw['snapObjs']).max() for p in prefxs]),
    })

def dump_summary(system : str = SYSTMS[0], sigfig : int = 6) -> None:
    """
    Saves the summarized results as a LaTeX-ready table
    """
    system = system.upper()
    assert system in SYSTMS

    raw = get_raw_result(system=system, drop_infeasible=True)
    tab = tabulate_result(raw)

    FNAME = f"summary_{system}"
    LCAPT = f"Results for {system}. There are {raw['areFeasible'].sum()} instances."
    SCAPT = f"Results for {system}."
    tab.to_latex(
        buf=f"./summary_{args.system.upper()}.tex",
        index=False, escape=False, caption=(LCAPT, SCAPT),
        formatters={0 : _float2sci},
        float_format=f"%.{sigfig}g",
    )
    print(f"{' ' * 3}Table saved to '{FNAME}.tex'")

if __name__ == "__main__":
    args = get_argparser().parse_args()

    print(f"Summarizing results for {args.system.upper()}")
    dump_summary(system=args.system, sigfig=args.sigfig)
