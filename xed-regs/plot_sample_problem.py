"""
A script plotting the illustrative problem

Masterarbeit-DemApp
Copyright (C) 2022 - present, Christian Cahig
https://github.com/christian-cahig/Masterarbeit-DemApp

This file is part of Masterarbeit-DemApp, and is covered by the CC-BY-4.0 License.
See `LICENSE` file in the source directory for details.
"""

from argparse import ArgumentParser, Namespace

import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
import seaborn as sns

__author__ = "Christian Cahig"

def get_argparser() -> ArgumentParser:
    """
    Initializes the argument parser
    """
    ap = ArgumentParser(description="A script for plotting the illustrative problem")
    ap.add_argument(
        "--xmin", type=float, required=False, default=-3.,
        help="Lower bound for the decision variable",
    )
    ap.add_argument(
        "--xmax", type=float, required=False, default=3.,
        help="Upper bound for the decision variable",
    )
    ap.add_argument(
        "--Nx", type=int, required=False, default=1_000,
        help="Number of samples drawn from the domain [XMIN, XMAX]",
    )
    ap.add_argument(
        "--rmin", type=float, required=False, default=0.,
        help="Minimum value for the regularization weight",
    )
    ap.add_argument(
        "--rmax", type=float, required=False, default=2.,
        help="Maximum value for the regularization weight",
    )
    ap.add_argument(
        "--Nr", type=int, required=False, default=5,
        help="Number of evenly-spaced values drawn from the interval [RMIN, RMAX]",
    )
    return ap

def objfun(
    dom : NDArray,
    reg : float = 0.,
) -> NDArray:
    """
    Evaluates the objective function given the variable value and regularization
    """
    return np.abs(dom) + (reg * np.abs(dom - 1))

def main(
    xlims : tuple[float, float],
    Nx : int = 1_000,
    rlims : tuple[float, float] = (0., 2.),
    Nr : int = 5,
) -> None:
    """
    Main routine
    """
    dom = np.linspace(min(xlims), max(xlims), num=Nx)
    regs = np.linspace(min(rlims), max(rlims), num=Nr)

    FNAME = "sample"
    fig = plt.figure(figsize=(6.01, 2.7), dpi=100)
    gsp = fig.add_gridspec(nrows=1, ncols=1)
    axs = fig.add_subplot(gsp[0, 0])

    with sns.axes_style("whitegrid"):
        for reg in regs:
            sns.lineplot(
                x=dom, y=objfun(dom, reg=reg), ax=axs,
                palette="pastel",
                linewidth=0.7,
                label=f"{reg:.1f}",
            )

    axs.grid(visible=True, axis="both", which="major")
    axs.set_xlim(left=min(xlims) - 0.1, right=max(xlims) + 0.1)
    axs.set_ylim(bottom=-0.1, top=3.1)
    axs.set_xlabel(r"Variable, $x$", fontsize=8)
    axs.set_ylabel(r"Objective, $f$", fontsize=8)
    for l in axs.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs.get_yaxis().get_ticklabels(): l.set_fontsize(7)
    axs.legend(
        loc="best", fontsize=7, ncol=2,
        title=r"Regularization, $\mu$",
        title_fontproperties={"size" : 8},
    )

    fig.tight_layout()
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"Plots saved to '{FNAME}.pdf'")

if __name__ == "__main__":
    args = get_argparser().parse_args()
    main(
        (args.xmin, args.xmax), Nx=args.Nx,
        rlims=(args.rmin, args.rmax), Nr=args.Nr,
    )

