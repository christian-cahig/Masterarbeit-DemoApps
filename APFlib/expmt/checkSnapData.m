function [flag, acct] = checkSnapData(sdata, NameVals)
%CHECKSNAPDATA Sanity checks for snapshot data
%   [flag, acct] = CHECKSNAPDATA(sdata, NameVals)
%   performs a series of sanity checks to determine the validity of a given snapshot data set.
%
%   Inputs
%   ---------------------
%   sdata : the MATPOWER case struct (in MATPOWER's internal indexing) having the snapshot data.
%   NameVals : optional name-value pairs.
%   (-) EqualityTol, a small positive value used as the "small enough" tolerance when checking for
%       equality between quantities (e.g., norm of two vectors expected to be equal). When unspeci-
%       fied, the tolerance used is 1e-3.
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
%   https://doi.org/10.5281/zenodo.7523544
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    sdata {mustBeMPC}
    NameVals.EqualityTol (1, 1) {mustBePositive} = 1e-3
    NameVals.beVerbose (1, 1) {mustBeNumericOrLogical} = true
end

if ~isInterIdxMPC(sdata), sdata = ext2int(sdata); end

%% Utilities
INDENT = "   ";
NORM_ORDERS = [2; Inf];
NORM_ORDERS_STR = {'euc'; 'inf'};
NORM_NOT_SMOL = "> eps";

%% Intro
if NameVals.beVerbose, fprintf("%sSanity-checking snapshot data.\n", INDENT); end
acct = {};

[Cf, Ct] = makeCfCt(sdata);
Cu = makeCu(sdata);
Cd = makeCd(sdata);
Y = makeY(sdata);
ref_bus = bustypes(sdata.bus, sdata.gen);
% are_nonref_buses = (1 : N_b) ~= ref_bus;
N_b = size(sdata.bus, 1);
% N_br = size(sdata.branch, 1);
% N_u = size(sdata.gen, 1);
[dem_buses, ~] = getDemBuses(sdata);
vm = sdata.bus(:, 8);
va = deg2rad(sdata.bus(:, 9));
pu = sdata.gen(:, 2);
qu = sdata.gen(:, 3);
pd = sdata.bus(dem_buses, 3);
qd = sdata.bus(dem_buses, 4);
clear dem_buses;

%% Feasibility of supply injections
if NameVals.beVerbose, fprintf("%s%sFeasibility of supply injections...", INDENT, INDENT); end
if any(pu > sdata.gen(:, 9))
    acct = [acct; {"pu > pu_max"}];
end
if any(pu < sdata.gen(:, 10))
    acct = [acct; {"pu < pu_min"}];
end
if any(qu > sdata.gen(:, 4))
    acct = [acct; {"qu > qu_max"}];
end
if any(qu < sdata.gen(:, 5))
    acct = [acct; {"qu < qu_min"}];
end
if NameVals.beVerbose, fprintf("done.\n"); end

%% Feasibility of bus voltage magnitudes
if NameVals.beVerbose, fprintf("%s%sFeasibility of bus voltage magnitudes...", INDENT, INDENT); end
if any(sdata.bus(:, 13) <= 0)
    acct = [acct; {"vm_min <= 0"}];
end
if any(sdata.bus(:, 12) <= 0)
    acct = [acct; {"vm_max <= 0"}];
end
if any(sdata.bus(:, 13) > sdata.bus(:, 12))
    acct = [acct; {"vm_min > vm_max"}];
end
if any(vm < sdata.bus(:, 13))
    acct = [acct; {"vm < vm_min"}];
end
if any(vm > sdata.bus(:, 12))
    acct = [acct; {"vm < vm_max"}];
end
if NameVals.beVerbose, fprintf("done.\n"); end

%% Net nodal injections
if NameVals.beVerbose, fprintf("%s%sNet nodal injections...", INDENT, INDENT); end
e_cd = calcPQFromCD(Cu, Cd, ...
    pu ./ sdata.gen(:, 7), qu ./ sdata.gen(:, 7), ...
    pd ./ sdata.baseMVA, qd ./ sdata.baseMVA ...
); 
clear Cu Cd pu qu pd qd;
% 
[e_sY, J] = calcPQJacPQ_VmVa(vm, va, Y);
clear vm va Y;
% 
for i = 1:length(NORM_ORDERS)
    if norm(e_sY - e_cd, NORM_ORDERS(i)) > NameVals.EqualityTol
        msg = sprintf("%snorm(e(s) - e(c,d)) %s", NORM_ORDERS_STR{i}, NORM_NOT_SMOL);
        acct = [acct; {msg}];
    end
end
if NameVals.beVerbose, fprintf("done.\n"); end

%% Nodal power residuals
if NameVals.beVerbose, fprintf("%s%sNodal power residuals...", INDENT, INDENT); end
phi = e_sY - e_cd;
clear e_cd e_sY;
% 
for i = 1:length(NORM_ORDERS)
    if norm(J' * phi, NORM_ORDERS(i)) > NameVals.EqualityTol
        msg = sprintf("%snorm(grad(0.5*eucnorm(phi(s,c,d))^2, s)) %s", NORM_ORDERS_STR{i}, NORM_NOT_SMOL);
        acct = [acct; {msg}];
    end
end
clear J phi;
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
