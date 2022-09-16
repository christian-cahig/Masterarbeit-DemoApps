"""
A script summarizing the run time metrics for calculating the anticipated bus voltages

Masterarbeit-DemApp
Copyright (C) 2022 - present, Christian Cahig
https://github.com/christian-cahig/Masterarbeit-DemApp

This file is part of Masterarbeit-DemApp, and is covered by the CC-BY-4.0 License.
See `LICENSE` file in the source directory for details.
"""

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
from collections import Counter

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.io import loadmat
import matplotlib.pyplot as plt
from numpy.typing import NDArray

__author__ = "Christian Cahig"
__version__ = "0.1.0"
__all__ = []

SYSTMS = ["AEP57", "PGE69", "PEG89", "EDC141", "ACT500", "RTE1888", "POL3375"]
SOLVRS = ["Powell Hybrid", "Levenberg-Marquardt"]
SOLVID = ["PH", "LM"]
TSCALE = 1e-3
TPREFX = "m"

def float_to_scinotation_str(
    num : float,
    digits : int = 3,
    format : str = "g",
) -> str:
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

def get_raw_result(
    mat_file : str,
    solver : str = SOLVRS[0],
    drop_infeasible : bool = False,
) -> dict[str, NDArray]:
    """
    Loads the raw results concerning one system
    """
    assert solver in SOLVRS

    out = loadmat(mat_file, squeeze_me=True, simplify_cells=True)
    out['SOLVER'] = solver
    solver = SOLVID[SOLVRS.index(solver)].lower()

    # Remove unnecessary entries
    for k in [k for k in out if (
        (not k.startswith(solver))
        and (
            k.endswith("Oks")
            or k.endswith("Norms")
            or k.endswith("Evals")
            or k.endswith("Steps")
            or k.endswith("Times")
            or k.endswith("__")
        )
    )]: del out[k]

    # Convert to booleans
    out['areFeasible'] = out['areFeasible'].astype(bool)
    out[f"{solver}Oks"] = out[f"{solver}Oks"].astype(bool)

    # Convert to integers
    for k in [k for k in out if (k.endswith("Steps") or k.endswith("Evals"))]:
        out[k] = out[k].astype(np.uint8)

    # Drop power-flow-infeasible instances
    if drop_infeasible:
        for k in [k for k in out if k.startswith(solver)]:
            out[k] = out[k][out['areFeasible']]
        out['areFeasible'] = out['areFeasible'][out['areFeasible']]

    # Rename
    for k in [k for k in out if k.startswith(solver)]:
        out[f"apf{k.split(solver)[-1]}"] = out[k]
        del out[k]

    return out

def prep_results(
    systems: list[str] = SYSTMS,
    solver: str = SOLVRS[0],
) -> pd.DataFrame:
    """
    Loads and preprocesses the results concerning all systems
    """
    assert solver in SOLVRS

    raw = {}
    for sys in sorted(systems):
        res = get_raw_result(f"{sys}_results.mat", solver=solver, drop_infeasible=True)
        raw[res['SNAPSHOT']] = res
    del sys, res

    num_solved = [raw[k]['apfOks'].sum() for k in raw]

    df = pd.DataFrame(data={
        "Iterations" : pd.Series(data=np.concatenate(
            [raw[k]['apfSteps'] for k in raw],
            axis=0,
        )),
        "Func. evals" : pd.Series(data=np.concatenate(
            [raw[k]['apfEvals'] for k in raw],
            axis=0
        )),
        "Time" : pd.Series(data=np.concatenate(
            [raw[k]['apfTimes'] / TSCALE for k in raw],
            axis=0,
        )),
        "System" : pd.Series(
            data=list(Counter(dict(zip(systems, num_solved))).elements()),
            dtype="string",
        ),
    })
    del raw

    return df

def summarize_results(df : pd.DataFrame) -> pd.DataFrame:
    """
    Summarizes, in terms of extrema, means, and standard deviations, the results for each system
    """
    out = df.groupby(["System"]).count().reindex(SYSTMS).rename_axis(None, axis="rows")
    out.drop(columns=["Func. evals", "Time"], inplace=True)
    out.rename(columns={out.columns[0] : "Instances"}, inplace=True)

    _cols = f"steps|evals|time [{TPREFX}s]".split("|")
    _mins = df.groupby(["System"]).min().reindex(SYSTMS)
    _mins.rename(columns={c : f"Min. {c_}" for c, c_ in zip(_mins.columns, _cols)}, inplace=True)
    _maxs = df.groupby(["System"]).max().reindex(SYSTMS)
    _maxs.rename(columns={c : f"Max. {c_}" for c, c_ in zip(_maxs.columns, _cols)}, inplace=True)
    _avgs = df.groupby(["System"]).mean().reindex(SYSTMS)
    _avgs.rename(columns={c : f"Avrg. {c_}" for c, c_ in zip(_avgs.columns, _cols)}, inplace=True)
    _sdvs = df.groupby(["System"]).std().reindex(SYSTMS)
    _sdvs.rename(columns={c : f"Stdv. {c_}" for c, c_ in zip(_sdvs.columns, _cols)}, inplace=True)

    out = pd.concat((
        out,
        _mins[_mins.columns[0]], _maxs[_maxs.columns[0]],
        _avgs[_avgs.columns[0]], _sdvs[_sdvs.columns[0]],
        _mins[_mins.columns[1]], _maxs[_maxs.columns[1]],
        _avgs[_avgs.columns[1]], _sdvs[_sdvs.columns[1]],
        _mins[_mins.columns[2]], _maxs[_maxs.columns[2]],
        _avgs[_avgs.columns[2]], _sdvs[_sdvs.columns[2]],
    ), axis=1)
    out.index.name = "System"

    return out.reset_index()

def summarize_solver_results(
    solver : str = SOLVRS[0],
    niter_xlim : tuple[int | float | None, int | float | None] = (None, None),
    neval_xlim : tuple[int | float | None, int | float | None] = (None, None),
    wtime_xlim : tuple[int | float | None, int | float | None] = (None, None),
) -> None:
    """
    Summarizes and plots the results for a specified solver
    """
    assert solver in SOLVRS
    solver_ = SOLVID[SOLVRS.index(solver)].upper()

    df = prep_results(systems=SYSTMS, solver=solver)
    print(f"Summarizing results for {solver}")

    # Summary table
    FNAME = f"apf2-summary_{solver_}"
    LONGCAPT = f"This is the long caption for {solver}."
    SHORCAPT = f"A short caption for {solver}."
    dfs = summarize_results(df)
    dfs.to_latex(
        buf=f"./{FNAME}.tex",
        index=False, escape=False, caption=(LONGCAPT, SHORCAPT),
        formatters={0 : float_to_scinotation_str},
    )
    print(f"{' ' * 3}Summary table saved to '{FNAME}.tex'")
    del dfs

    # Solver iterations
    FNAME = f"apf2-niters_{solver_}"
    fig_steps = plt.figure(figsize=(6.0, 3.5), dpi=100)
    gsp_steps = fig_steps.add_gridspec(nrows=1, ncols=1)
    axs_steps = fig_steps.add_subplot(gsp_steps[0, 0])
    with sns.axes_style("whitegrid"):
        axs_steps = sns.stripplot(
            x="Iterations", y="System", data=df, order=SYSTMS,
            marker=".", size=10, jitter=0.3,
            palette="pastel", linewidth=0.1,
        )
    axs_steps.set_xlim(left=niter_xlim[0], right=niter_xlim[1])
    axs_steps.grid(visible=True, axis="x", which="both")
    axs_steps.set_xlabel(
        f"Number of {solver} iterations",
        fontsize=8,
    )
    axs_steps.set_ylabel("")
    for l in axs_steps.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_steps.get_yaxis().get_ticklabels(): l.set_fontsize(8); l.set_fontfamily("monospace")
    fig_steps.tight_layout()
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots of solver iterations saved to '{FNAME}.pdf'")
    del fig_steps, gsp_steps, axs_steps

    # Solver function evaluations
    FNAME = f"apf2-nevals_{solver_}"
    fig_evals = plt.figure(figsize=(6.0, 3.5), dpi=100)
    gsp_evals = fig_evals.add_gridspec(nrows=1, ncols=1)
    axs_evals = fig_evals.add_subplot(gsp_evals[0, 0])
    with sns.axes_style("whitegrid"):
        axs_evals = sns.stripplot(
            x="Func. evals", y="System", data=df, order=SYSTMS,
            marker=".", size=10, jitter=0.3,
            palette="pastel", linewidth=0.1,
        )
    axs_evals.set_xlim(left=neval_xlim[0], right=neval_xlim[1])
    axs_evals.grid(visible=True, axis="x", which="both")
    axs_evals.set_xlabel(
        f"Number of merit function evaluations via {solver}",
        fontsize=8,
    )
    axs_evals.set_ylabel("")
    for l in axs_evals.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_evals.get_yaxis().get_ticklabels(): l.set_fontsize(8); l.set_fontfamily("monospace")
    fig_evals.tight_layout()
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots of solver merit function evaluations saved to '{FNAME}.pdf'")
    del fig_evals, gsp_evals, axs_evals

    # Solver run times
    FNAME = f"apf2-wtimes_{solver_}"
    fig_time = plt.figure(figsize=(6.0, 3.5), dpi=100)
    gsp_time = fig_time.add_gridspec(nrows=1, ncols=1)
    axs_time = fig_time.add_subplot(gsp_time[0, 0])
    with sns.axes_style("whitegrid"):
        axs_time = sns.stripplot(
            x="Time", y="System", data=df, order=SYSTMS,
            marker=".", size=10, jitter=0.3,
            palette="pastel", linewidth=0.1,
        )
    axs_time.set_xlim(left=wtime_xlim[0], right=wtime_xlim[1])
    axs_time.grid(visible=True, axis="x", which="both")
    axs_time.set_xlabel(
        f"{solver} run time [{TPREFX}s]",
        fontsize=8,
    )
    axs_time.set_ylabel("")
    for l in axs_time.get_xaxis().get_ticklabels(): l.set_fontsize(7)
    for l in axs_time.get_yaxis().get_ticklabels(): l.set_fontsize(8); l.set_fontfamily("monospace")
    fig_time.tight_layout()
    plt.savefig(f"./{FNAME}.pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    print(f"{' ' * 3}Plots of wall times saved to '{FNAME}.pdf'")
    del fig_time, gsp_time, axs_time

if __name__ == "__main__":
    niter_xlims = [(0, 40), (0, 40)]
    neval_xlims = [(0, 60), (0, 60)]
    wtime_xlims = [(0, 400), (0, 400)]
    for s, n, f, t  in zip(SOLVRS, niter_xlims, neval_xlims, wtime_xlims):
        summarize_solver_results(solver=s, niter_xlim=n, neval_xlim=f, wtime_xlim=t)
