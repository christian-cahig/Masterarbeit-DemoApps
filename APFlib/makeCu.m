function [Cu] = makeCu(mpc)
%MAKECG Builds the supply unit connection matrix
%   Cu = MAKECG(mpc)
%   returns the supply unit connection matrix from a given MATPOWER case struct, assumed to be
%   using internal indexing.
%
%   With N_b and N_u as the respective the numbers of buses and of supply units, Cu is sparse with
%   N_b rows and N_u columns, where Cu(i,j) is 1 if the j-th supply unit is at the i-th bus.
%
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format and indexing, please refer to Appendix B and Section 9.4.
%   2. The main logic herein is based on MATPOWER's `bustypes` function.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7180586
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeInterIdxMPC}
end

%% Main
N_b = size(mpc.bus, 1);   % Number of buses
N_u = size(mpc.gen, 1);   % Number of supply units
Cu = sparse(mpc.gen(:,1), (1:N_u)', mpc.gen(:,8) > 0, N_b, N_u);

end
