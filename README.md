# *Demonstrations and Applications*

This repository houses the code files for
Chapter 4 (*Demonstrations and Applications*)
of my master's thesis
(*A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power Systems Analysis*).
The complete repository
(which includes the manuscript and the data files)
will be made available at
[`https://doi.org/10.5281/zenodo.7077324`](https://doi.org/10.5281/zenodo.7077324).

<!-- omit in toc -->
## Contents

- [Overview](#overview)
- [Citing](#citing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Overview

The core anticipatory power flow (APF) routines are packaged as a small MATLAB library [`APFlib`](./APFlib/).
This library is built on [MATPOWER](https://github.com/MATPOWER/matpower)
(specifically, [version 7.1](https://github.com/MATPOWER/matpower/releases/tag/7.1))
and uses [CVX](http://cvxr.com/cvx)
(specifically, [version 2.2 build 1148 (62bfcca)](http://cvxr.com/cvx/download/)).
For a complete list of the dependencies, please see Table 4.1 of the manuscript.
To "install" `APFlib`,

1. make sure that MATPOWER version 7.1 and CVX version 2.2 are available;
   and
2. with the current working directory at the root of this repository,
   invoke `installAPFlib` in the MATLAB Command Window to install `APFlib`.

Code files are organized as follows.

- Those pertaining to Section 4.1.1 (*Solving the extended economic dispatch subproblem*)
  are in `time-apf1/`.
- Those pertaining to Section 4.1.2 (*Solving the APF equations*) are in `time-apf2/`.
- Those pertaining to Section 4.2 (*Effects of the supply regularization strengths*)
  are in `xed-regs/`.
- Those pertaining to Section 4.3.1 (*Providing warm-start points*) are in `warm-start/`.
- Those pertaining to Section 4.3.2 (*Finding the nearest power-flow feasible point*)
  are in `nearest-feasible/`.
- Those pertaining to Section 4.4 (*Differentiating through the APF equations*)
  are in `diff-apf2/`.
- Those pertaining to Appendix B.2 (*The snapshot data sets*) are in `snapshots/`.

All experiments and examples are implemented as MATLAB scripts (*i.e.*, `*.m` files).
Snapshot data and raw results are available as `.mat` files
but **are not** included in this repository due to their size;
they can be accessed through the complete repository at
[`https://doi.org/10.5281/zenodo.7077324`](https://doi.org/10.5281/zenodo.7077324).

Processing of results (*e.g.*, producing the plots used in the manuscript)
are implemented as Python scripts.
To run these processing scripts,
instantiate a [Conda](https://docs.conda.io/projects/conda/en/latest/) environment
using the file [`env.yml`](./env.yml).
Please refer to the [Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html) for more details.

## Citing

We request that publications based on or derived from this repository,
or any part of it (*e.g.*, the `APFlib` library) or of my thesis,
explicitly acknowledge that fact by citing my work.
Please use (also available in [`CITATION.bib`](./CITATION.bib)):

```bibtex
@MastersThesis{CahigMastersThesis,
  author      = {Christian Y. {Cahig}},
  date        = {2022},
  institution = {Mindanao State University - Iligan Institute of Technology},
  title       = {A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power Systems Analysis},
  doi         = {10.5281/zenodo.7077324},
  address     = {Iligan, Philippines},
}
```

## License

This repository is licensed under the [Creative Commons Attribution 4.0 International Public License](https://creativecommons.org/licenses/by/4.0/).
Please see [`LICENSE`](./LICENSE) for the details.

## Acknowledgements

This work is supervised jointly by Dr. Marven E. Jabian and Dr. Abdul Aziz G. Mabaning,
and has been supported by the DOST-SEI ERDT program.
