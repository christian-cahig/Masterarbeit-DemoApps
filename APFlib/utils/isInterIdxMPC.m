function [flag] = isInterIdxMPC(mpc)
%ISINTERIDXMPC Checks if a MATPOWER case is internally numbered
%   flag = ISINTERIDXMPC(mpc)
%   returns a logical, indicating if `mpc` is an internally-ordered MATPOWER case.
%
%   There are five conditions that must be met for `mpc` to be considered internally indexed.
%      (1) There are no isolated buses in the bus data matrix.
%      (2) All branches are in service.
%      (3) All supply units online.
%      (4) The bus indices are numbered consecutively, starting at 1.
%      (5) The case struct `mpc` has a field `order`.
%
%   Notes
%   ---------------------
%   This function relies heavily on the MATPOWER package. It has been developed and tested on
%   version 7.1. The accompanying User's Manual is freely available at
%   https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format and on its indexing schemes, please refer to Appendix B and Section 9.4.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7491508
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeMPC}
end

%% Main
% Condition 1
if any(mpc.bus(:, 2) == 4)
    flag = false;
    return
end

% Condition 2
if any(mpc.branch(:, 11) < 1)
    flag = false;
    return
end

% Condition 3
if any(mpc.gen(:, 8) < 1)
    flag = false;
    return
end

% Condition 4
if any(mpc.bus(:, 1) ~= (1 : size(mpc.bus, 1))')
    flag = false;
    return
end

% Condition 5
if ~isfield(mpc, 'order')
    flag = false;
    return
end

flag = true;

end
