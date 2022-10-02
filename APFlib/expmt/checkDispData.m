function [flag, acct] = checkDispData(ddata, NameVals)
%CHECKDISPDATA Sanity checks for dispatch data
%   [flag, acct] = CHECKDISPDATA(ddata, NameVals)
%   performs a series of sanity checks to determine the validity of a given dispatch data set.
%
%   Inputs
%   ---------------------
%   ddata : the MATPOWER case struct (in MATPOWER's internal indexing) having the dispatch data.
%   NameVals : optional name-value pairs.
%   (-) beVerbose, a logical flag (`true` by default) for enabling verbosity.
%
%   Outputs
%   ---------------------
%   flag : a logical, indicating if the given dispatch data set passes all sanity checks.
%   acct : an optional cell array of strings, each indicating a test not passed.
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
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    ddata {mustBeMPC}
    NameVals.beVerbose (1, 1) {mustBeNumericOrLogical} = true
end

if ~isInterIdxMPC(ddata), ddata = ext2int(ddata); end

%% Utilities
INDENT = "   ";

%% Intro
if NameVals.beVerbose, fprintf("%sSanity-checking dispatch data.\n", INDENT); end
acct = {};

% [Cf, Ct] = makeCfCt(ddata);
% Cu = makeCu(ddata);
% Cd = makeCd(ddata);
% Y = makeY(ddata);
% ref_bus = bustypes(ddata.bus, ddata.gen);
% are_nonref_buses = (1 : N_b) ~= ref_bus;
% N_b = size(ddata.bus, 1);
% N_br = size(ddata.branch, 1);
% N_u = size(ddata.gen, 1);
[dem_buses, ~] = getDemBuses(ddata);
% vm = ddata.bus(:, 8);
% va = deg2rad(ddata.bus(:, 9));
% pu = ddata.gen(:, 2);
% qu = ddata.gen(:, 3);
pd = ddata.bus(dem_buses, 3);
qd = ddata.bus(dem_buses, 4);

%% Scheduled supply injection limits
if NameVals.beVerbose, fprintf("%s%sGenerator injection bounds...", INDENT, INDENT); end
pu_min = ddata.gen(:, 10);
pu_max = ddata.gen(:, 9);
qu_min = ddata.gen(:, 5);
qu_max = ddata.gen(:, 4);
if all(pu_max == 0)
    acct = [acct; {"pu_max = 0"}];
end
if any(pu_min > pu_max)
    acct = [acct; {"pu_min > pu_max"}];
end
if all(pu_min == pu_max)
    acct = [acct; {"pu_min = pu_max"}];
end
if any(qu_min > qu_max)
    acct = [acct; {"qu_min > qu_max"}];
end
if all(qu_min == qu_max)
    acct = [acct; {"qu_min = qu_max"}];
end
clear pu_min pu_max qu_min qu_max;
if NameVals.beVerbose, fprintf("done.\n"); end

%% Bus voltage magnitude bounds
if NameVals.beVerbose, fprintf("%s%sBus voltage magnitude bounds...", INDENT, INDENT); end
vm_min = ddata.bus(:, 13);
vm_max = ddata.bus(:, 12);
if any(vm_min <= 0)
    acct = [acct; {"vm_min <= 0"}];
end
if any(vm_max <= 0)
    acct = [acct; {"vm_max <= 0"}];
end
if any(vm_min > vm_max)
    acct = [acct; {"vm_min > vm_max"}];
end
clear vm_min vm_max;
if NameVals.beVerbose, fprintf("done.\n"); end

%% Anticipated load draws
if NameVals.beVerbose, fprintf("%s%sAnticipated load draws...", INDENT, INDENT); end
if any(pd < 0)
    acct = [acct; {"pd < 0"}];
end
clear pd qd;
if NameVals.beVerbose, fprintf("done.\n"); end

%% Summary
num_failed_checks = size(acct, 1);
if num_failed_checks == 0
    flag = true;
    summary_msg = sprintf("All checks passed :D");
else
    flag = false;
    summary_msg = sprintf("Failed %d check/s :(", num_failed_checks);
end

if NameVals.beVerbose, fprintf("%s%s\n", INDENT, summary_msg); end

end
