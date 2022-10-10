function [mpopt] = getGoToMPOpts()
%GETGOTOMPOPTS Loads a set of go-to MATPOWER options
%   mpopt = GETGOTOMPOPTS()
%   returns a MATPOWER options struct.
%
%   Notes
%   ---------------------
%   This function relies heavily on the MATPOWER package. It has been developed and tested on
%   version 7.1. The accompanying User's Manual is freely available at
%   https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details regarding MATPOWER options,
%   please refer to Appendix C.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7180586
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Main
mpopt = mpoption;

%% Verbosity and pretty-printing
mpopt.verbose = 0;
mpopt.out.all = 0;

%% PF
mpopt.pf.alg = 'NR';
mpopt.pf.nr.max_it = 100;

%% OPF formulation
mpopt.opf.flow_lim = 'S';
mpopt.opf.ignore_angle_lim = 1;

%% OPF solver
mpopt.opf.ac.solver = 'MIPS';
mpopt.opf.start = 2;
mpopt.mips.linsolver = '';
mpopt.mips.max_it = 300;
mpopt.mips.step_control = 1;

end
