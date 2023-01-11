function [Cd] = makeCd(mpc)
%MAKECD Builds the demand unit connection matrix.
%   Cd = MAKECD(mpc)
%   returns the demand unit connection matrix from a given MATPOWER case struct, assumed to be
%   using internal indexing.
%
%   With N_b and N_d as the respective numbers of buses and of demand units, Cd is a sparse matrix
%   with N_b rows and N_d matrix, with Cd(i,j) is 1 if the j-th demand unit is at the i-th bus.
%
%   Notes
%   ---------------------
%   This function relies heavily on the MATPOWER package. It has been developed and tested on
%   version 7.1. The accompanying User's Manual is freely available at
%   https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format and indexing, please refer to Appendix B and Section 9.4.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7523544
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeInterIdxMPC}
end

%% Main
N_b = size(mpc.bus, 1);                 % Number of buses
[LD_buses, N_d] = getDemBuses(mpc);  % Demand buses and number of demand units
Cd = sparse(LD_buses, (1:N_d)', true(N_d, 1), N_b, N_d);

end
