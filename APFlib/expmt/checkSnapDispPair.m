function [flag, acct] = checkSnapDispPair(sdata, ddata, NameVals)
%CHECKSNAPDISPPAIR Sanity checks for a pair of snapshot and dispatch data sets
%   [flag, acct] = CHECKSNAPDISPPAIR(sdata, ddata, NameVals)
%   performs a series of sanity checks to determine if the given sets of snapshot and dispatch
%   data are valid.
%
%   Inputs
%   ---------------------
%   sdata : the MATPOWER case struct (in MATPOWER's internal indexing) having the snapshot data.
%   ddata : the MATPOWER case struct (in MATPOWER's internal indexing) having the dispatch data.
%   NameVals : optional name-value pairs.
%   (-) beVerbose, a logical flag (`true` by default) for enabling verbosity.
%
%   Outputs
%   ---------------------
%   flag : a logical, indicating if the given snapshot data set passes all sanity checks.
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
%   https://doi.org/10.5281/zenodo.7180586
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    sdata {mustBeMPC}
    ddata {mustBeMPC}
    NameVals.beVerbose (1, 1) {mustBeNumericOrLogical} = true
end

if ~isInterIdxMPC(sdata), sdata = ext2int(sdata); end
if ~isInterIdxMPC(ddata), ddata = ext2int(ddata); end

%% Utilities
INDENT = "   ";

%% Intro
if NameVals.beVerbose, fprintf("%sSanity-checking snapshot and dispatch data.\n", INDENT); end
acct = {};

%% Base MVAs
if NameVals.beVerbose, fprintf("%s%sChecking base MVAs...", INDENT, INDENT); end
if sdata.baseMVA ~= ddata.baseMVA
    acct = [acct; {"Sys. base MVAs don't match"}];
end
if any(sdata.gen(:, 7) ~= ddata.gen(:, 7))
    acct = [acct; {"Sup. base MVAs don't match"}];
end
if NameVals.beVerbose, fprintf("done.\n"); end

%% Load units
if NameVals.beVerbose, fprintf("%s%sChecking demand buses...", INDENT, INDENT); end
snap_dem_buses = getDemBuses(sdata);
snap_N_d = numel(snap_dem_buses);
disp_dem_buses = getDemBuses(ddata);
disp_N_d = numel(disp_dem_buses);
if (snap_N_d ~= disp_N_d) || any(snap_dem_buses ~= disp_dem_buses)
    acct = [acct; {"Loaded buses don't match"}];
end
clear snap_dem_buses disp_dem_buses;
if NameVals.beVerbose, fprintf("done.\n"); end

%% Element counts
if NameVals.beVerbose, fprintf("%s%sChecking element counts...", INDENT, INDENT); end
if size(sdata.bus, 1) ~= size(ddata.bus, 1)
    acct = [acct; {"Bus counts don't match"}];
end
if size(sdata.branch, 1) ~= size(ddata.branch, 1)
    acct = [acct; {"Branch counts don't match"}];
end
if size(sdata.gen, 1) ~= size(ddata.gen, 1)
    acct = [acct; {"Sup. unit counts don't match"}];
end
if snap_N_d ~= disp_N_d
    acct = [acct; {"Load unit counts don't match"}];
end
if NameVals.beVerbose, fprintf("done.\n"); end

%% Summary
num_failed_checks = size(acct, 1);
if num_failed_checks > 0
    flag = false;
    summary_msg = sprintf("Failed %d check/s :(", num_failed_checks);
else
    flag = true;
    summary_msg = sprintf("All checks passed :D");
end

if NameVals.beVerbose, fprintf("%s%s\n", INDENT, summary_msg); end

end
