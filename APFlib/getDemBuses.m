function [dem_buses, N_d] = getDemBuses(mpc)
%GETDEMBUSES Returns the indices of demand buses
%   [dem_buses, N_d] = GETDEMBUSES(mpc)
%   returns the indices of buses with demand units, and optionally, the number of demand units,
%   from a given MATPOWER case struct, assumed to be using internal indexing.
%
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%      format and indexing, please refer to Appendix B and Section 9.4.
%   2. A 'demand bus' is defined here as that whose corresponding entries in columns 3 and 4 of the
%      bus data matrix, `mpc.bus`, are both nonzero.
%   3. Since MATPOWER does not have a dedicated data matrix for demand units, a demand bus is
%      assumed to cater one demand unit. Hence, the number of demand buses is taken to be the
%      number of demand units as well.
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
dem_buses = find((mpc.bus(:, 3) ~= 0) | (mpc.bus(:, 4) ~= 0));

if nargout > 1  % Number of demand units
    N_d = length(dem_buses);
end

end
