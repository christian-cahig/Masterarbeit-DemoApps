function [Cf, Ct] = makeCfCt(mpc)
%MAKECFCT Builds the branch connection matrices.
%   [Cf, Ct] = MAKECFCT(mpc)
%   returns the branch connection matrices from a given MATPOWER case struct, assumed to be using
%   internal indexing.
%
%   With N_b and N_br as the respective numbers of buses and of branches, Cf and Ct are sparse
%   N_br-by-N_b matrices, where Cf(i,j) and Ct(i,k) are 1 if the i-th branch connects the j-th bus
%   to the k-th bus.
%
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format and indexing, please refer to Appendix B and Section 9.4.
%   2. The main logic herein is based on MATPOWER's `makeYbus` function.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeInterIdxMPC}
end

%% Main
N_b = size(mpc.bus, 1);     % Number of buses
N_br = size(mpc.branch, 1); % Number of branches
f_buses = mpc.branch(:, 1); % "From"-side buses
t_buses = mpc.branch(:, 2); % "To"-side buses

Cf = sparse(1:N_br, f_buses, true(N_br, 1), N_br, N_b);
Ct = sparse(1:N_br, t_buses, true(N_br, 1), N_br, N_b);

end
