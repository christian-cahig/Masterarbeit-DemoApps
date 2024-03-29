function [lin, num] = getLines(mpc, NameVals)
%GETLINES Returns the indices of transmission lines or cables
%   [lin, num] = GETLINES(mpc)
%   returns the indices of branches representing transmission lines or cables, and optionally, the
%   number of them, given a MATPOWER case struct assumed to be using internal indexing.
%   [lin, num] = GETLINES(mpc, 'getFlowLtd', true)
%   returns the indices of branches representing transmission lines or cables with constraints on
%   power flows, and optionally, the number of them.
%
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%      format and indexing, please refer to Appendix B and Section 9.4.
%   2. 'Transmission lines' and 'cables' are defined here as branches whose corresponding column-9
%      entries in the branch data matrix, `mpc.branch`, are zero.
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
    NameVals.getFlowLtd (1, 1) {mustBeNumericOrLogical} = false
end

%% Main
if NameVals.getFlowLtd
    lin = find((mpc.branch(:, 9) == 0) & (mpc.branch(:, 6) > 0));
else
    lin = find(mpc.branch(:, 9) == 0);
end

if nargout > 1, num = length(lin); end

end
