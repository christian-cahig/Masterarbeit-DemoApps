function [flag, acct] = checkProbData(pdata, NameVals)
%CHECKPROBDATA Sanity checks for a given APF problem data
%   [flag, acct] = CHECKPROBDATA(problem, NameVals)
%   performs a series of sanity checks to determine the validity of a set of APF problem data.
%
%   Inputs
%   ---------------------
%   pdata : the problem data struct.
%   NameVals : optional name-value pairs.
%   (-) beVerbose, a logical flag (`true` by default) for enabling verbosity.
%
%   Outputs
%   ---------------------
%   flag : a logical, indicating if the given problem data set passes all sanity checks.
%   acct : an optional cell array of strings, each indicating a test not passed.
%
%   Notes
%   ---------------------
%   This function relies heavily on the MATPOWER package. It has been developed and tested on
%   version 7.1. The accompanying User's Manual is freely available at
%   https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format, please refer to Appendix B.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7523544
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    pdata
    NameVals.beVerbose (1, 1) {mustBeNumericOrLogical} = true
end

%% Utilities
INDENT = "   ";

%% Intro
if NameVals.beVerbose, fprintf("%sSanity-checking APF problem data.\n", INDENT); end
acct = {};

%% Snapshot data
if NameVals.beVerbose, fprintf("%s%sChecking snapshot data...", INDENT, INDENT); end
if ~isfield(pdata, 'snapshot')
    acct = [acct; {"Missing field `snapshot`"}];
elseif ~isMPC(pdata.snapshot)
    acct = [acct; {"Field `snapshot` is not a MATPOWER case"}];
elseif ~checkSnapData(pdata.snapshot, 'beVerbose', false)
    acct = [acct; {"Invalid snapshot data"}];
end
if NameVals.beVerbose, fprintf("done.\n"); end

%% Dispatch data
if NameVals.beVerbose, fprintf("%s%sChecking dispatch data...", INDENT, INDENT); end
if ~isfield(pdata, 'dispatch')
    acct = [acct; {"Missing field `dispatch`"}];
elseif ~isMPC(pdata.dispatch)
    acct = [acct; {"Field `dispatch` is not a MATPOWER case"}];
elseif ~checkDispData(pdata.dispatch, 'beVerbose', false)
    acct = [acct; {"Invalid dispatch data"}];
end
for f={'ref_bus', 'Cf', 'Ct', 'Y', 'va_ref'}
    if ~isfield(pdata, f{1})
        acct = [acct; {sprintf("Missing field '%s'", f{1})}];
    end
end
if NameVals.beVerbose, fprintf("done.\n"); end

%% Snapshot-dispatch pair
if NameVals.beVerbose, fprintf("%s%sChecking dispatch data...", INDENT, INDENT); end
if all(isfield(pdata, {'snapshot', 'dispatch'})) && ~checkSnapDispPair(pdata.snapshot, pdata.dispatch, 'beVerbose', false)
    acct = [acct; {"Invalid snapshot-dispatch pair"}];
end
if NameVals.beVerbose, fprintf("done.\n"); end

%% Invariant data
if NameVals.beVerbose, fprintf("%s%sChecking invariant data...", INDENT, INDENT); end
for f={'sys_baseMVA', 'sup_baseMVA', 'dem_buses', 'N_b', 'N_br', 'N_u', 'N_d', 'Cd', 'Cu'}
    if ~isfield(pdata, f{1})
        acct = [acct; {sprintf("Missing field '%s'", f{1})}];
    end
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
