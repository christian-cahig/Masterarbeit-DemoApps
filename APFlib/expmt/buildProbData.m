function [pdata] = buildProbData(sdata, ddata, NameVals)
%BUILDPROBDATA Sets up the problem data for an APF instance
%   pdata = BUILDPROBDATA(sdata, ddata, NameVals)
%   collects the APF problem data into one struct.
%
%   Inputs
%   ---------------------
%   sdata : the MATPOWER case struct (in MATPOWER's internal indexing) having the snapshot data.
%   ddata : the MATPOWER case struct (in MATPOWER's internal indexing) having the dispatch data.
%   NameVals : optional name-value pairs.
%   (-) calcCft, a logical flag (`false` by default), indicating whether or not the branch connec-
%       tion matrices are computed.
%   (-) calcCu, a logical flag (`false` by default), indicating whether or not the supply-unit con-
%       nection matrix is computed.
%   (-) calcCd, a logical flag (`false` by default), indicating whether or not the demand-unit con-
%       nection matrix is computed.
%   (-) calcY, a logical flag (`false` by default), indicating whether or not the system bus admit-
%       tance matrix is computed from the dispatch data.
%   (-) RefAngle, the reference value (in radians) for bus voltage phase angles. If unspecified, it
%       is set to 0.
%
%   Outputs
%   ---------------------
%   pdata : a struct containing the APF problem data. Its fields are as follows.
%   (-) snapshot, the snapshot data struct.
%   (-) dispatch, the dispatch data struct.
%   (-) sys_baseMVA, the value of `{sdata,ddata}.baseMVA`.
%   (-) sup_baseMVA, the value of `{sdata,ddata}.gen(:, 7)`.
%   (-) ref_bus, the index of the reference bus, set to be the slack bus in `ddata`.
%   (-) dem_buses, a vector of indices of demand buses in `ddata`.
%   (-) N_{b,br,g,d}, the number of {buses, branches, supply units, demand units}.
%   (-) C{u,d,f,t}, the {supply, load, from-side branch, to-side branch} connection matrix, set to
%       `NaN` if name-value argument `calcC{u,d,ft,ft}` is `false`.
%   (-) Y, the dispatch-data system bus admittance matrix, set to `NaN` if name-value argument
%       `calcY` is `false`.
%   (-) va_ref, the reference (in radians) for bus voltage phase angles.
%
%   Notes
%   ---------------------
%   This function relies heavily on the MATPOWER package. It has been developed and tested on
%   version 7.1. The accompanying User's Manual is freely available at
%   https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format and indexing, please refer to Appendix B and Section 9.4.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    sdata {mustBeMPC}
    ddata {mustBeMPC, mustBeValidSnapDispPair(sdata, ddata)}
    NameVals.calcCft (1, 1) {mustBeNumericOrLogical} = false
    NameVals.calcCu (1, 1) {mustBeNumericOrLogical} = false
    NameVals.calcCd (1, 1) {mustBeNumericOrLogical} = false
    NameVals.calcY (1, 1) {mustBeNumericOrLogical} = false
    NameVals.RefAngle (1, 1) {mustBeNumeric} = 0
end

%% Utilities
INDENT = "   ";

%% Intro
fprintf("%sBuilding problem data.\n", INDENT);
fprintf("%s%sCalc. Cf & Ct? %s\n", INDENT, INDENT, string(~~NameVals.calcCft));
fprintf("%s%sCalc. Cu? %s\n", INDENT, INDENT, string(~~NameVals.calcCu));
fprintf("%s%sCalc. Cd? %s\n", INDENT, INDENT, string(~~NameVals.calcCd));
fprintf("%s%sCalc. Y? %s\n", INDENT, INDENT, string(~~NameVals.calcY));
fprintf("%s%sRef. angle: %.3f rad.\n", INDENT, INDENT, NameVals.RefAngle);

% return
t_start = tic;

%% Snapshot data
fprintf("%s%sParsing snapshot data...", INDENT, INDENT);
if ~isInterIdxMPC(sdata)
    pdata.snapshot = ext2int(sdata);
else
    pdata.snapshot = sdata;
end
fprintf("done.\n");

%% Dispatch data
fprintf("%s%sParsing dispatch data...", INDENT, INDENT);
if ~isInterIdxMPC(ddata)
    pdata.dispatch = ext2int(ddata);
else
    pdata.dispatch = ddata;
end
pdata.ref_bus = bustypes(pdata.dispatch.bus, pdata.dispatch.gen);
pdata.are_nonref_buses = (1 : size(pdata.snapshot.bus, 1)) ~= pdata.ref_bus;
pdata.va_ref = NameVals.RefAngle;

if NameVals.calcCft
    [pdata.Cf, pdata.Ct] = makeCfCt(pdata.dispatch);
else
    [pdata.Cf, pdata.Ct] = deal(NaN, NaN);
end
if NameVals.calcY, pdata.Y = makeY(pdata.dispatch); else, pdata.Y = NaN; end

fprintf("done.\n");

%% Invariant data
fprintf("%s%sParsing invariant data...", INDENT, INDENT);
pdata.sys_baseMVA = pdata.dispatch.baseMVA;
pdata.sup_baseMVA = pdata.dispatch.gen(:, 7);
pdata.N_b = size(pdata.snapshot.bus, 1);
pdata.N_br = size(pdata.dispatch.branch, 1);
pdata.N_u = size(pdata.snapshot.gen, 1);
[pdata.dem_buses, pdata.N_d] = getDemBuses(pdata.dispatch);

if NameVals.calcCd, pdata.Cd = makeCd(pdata.dispatch); else, pdata.Cd = NaN; end
if NameVals.calcCu, pdata.Cu = makeCu(pdata.snapshot); else, pdata.Cu = NaN; end
fprintf("done.\n");

%% Summary
fprintf("%sProblem data built in %.5f s.\n", INDENT, toc(t_start));

end
