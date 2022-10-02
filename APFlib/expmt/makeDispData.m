function [ddata] = makeDispData(sdata, NameVals)
%MAKEDISPDATA Derives dispatch data from snapshot data
%   ddata = MAKEDISPDATA(sdata, NameVals)
%   synthesizes a set of dispatch data from a given snapshot data.
%
%   Inputs
%   ---------------------
%   sdata : a MATPOWER case struct (in MATPOWER's internal indexing) containing the snapshot data.
%   NameVals : optional name-value pairs.
%   (-) MPOptions, a struct containing the MATPOWER options. If unspecified, the one provided by
%       `getGoToMPOpts` is used.
%   (-) RScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.95; 1.05]` by default) for the branch series resistances in `sdata`.
%   (-) XScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.95; 1.05]` by default) for the branch series reactances in `sdata`.
%   (-) BScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.95; 1.05]` by default) for the branch shunt reactances in `sdata`.
%   (-) PdScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.80; 1.20]` by default) for the active demand draws in `sdata`.
%   (-) QdScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.80; 1.20]` by default) for the reactive demand draws in `sdata`.
%   (-) Seed, a nonnegative integer with which the random number generator is seeded. If unspecified
%       the random number generator is not seeded.
%   (-) assertOPFOk, a logical flag (`true` by default), indicating whether to require that solving
%       an OPF instance on the dispatch data succeeds.
%
%   Output
%   ---------------------
%   ddata : a MATPOWER case struct (in MATPOWER's internal indexing) containing the dispatch data.
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
    sdata {mustBeMPC}
    NameVals.MPOptions struct = getGoToMPOpts()
    NameVals.RScaleRange (2, 1) {mustBePositive} = [0.95; 1.05]
    NameVals.XScaleRange (2, 1) {mustBePositive} = [0.95; 1.05]
    NameVals.BScaleRange (2, 1) {mustBePositive} = [0.95; 1.05]
    NameVals.PdScaleRange (2, 1) {mustBePositive} = [0.80; 1.20]
    NameVals.QdScaleRange (2, 1) {mustBePositive} = [0.80; 1.20]
    NameVals.Seed {mustBeInteger, mustBeNonnegative, mustBeScalarOrEmpty} = []
    NameVals.assertOPFOk (1, 1) {mustBeNumericOrLogical} = true
end

%% Utilities
INDENT = "   ";
TRY_AGAIN_STR = "Try running again (if unseeded) or use a different seed.";
NO_ACOPF_SOLN_STR = sprintf("%s%sNo ACOPF solution found.\n%s%s\n", INDENT, INDENT, INDENT, TRY_AGAIN_STR);
if isempty(NameVals.Seed), SEED_STR = "None"; else, SEED_STR = sprintf("%d", NameVals.Seed); end

%% Intro
fprintf("%sBuilding dispatch data from a snapshot data set.\n", INDENT);
fprintf("%s%sBranch R scale: (%.3f, %.3f)\n", INDENT, INDENT, ...
                                              min(NameVals.RScaleRange), ...
                                              max(NameVals.RScaleRange));
fprintf("%s%sBranch X scale: (%.3f, %.3f)\n", INDENT, INDENT, ...
                                              min(NameVals.XScaleRange), ...
                                              max(NameVals.XScaleRange));
fprintf("%s%sBranch B scale: (%.3f, %.3f)\n", INDENT, INDENT, ...
                                              min(NameVals.BScaleRange), ...
                                              max(NameVals.BScaleRange));
fprintf("%s%sDemand P scale: (%.3f, %.3f)\n", INDENT, INDENT, ...
                                              min(NameVals.PdScaleRange), ...
                                              max(NameVals.PdScaleRange));
fprintf("%s%sDemand Q scale: (%.3f, %.3f)\n", INDENT, INDENT, ...
                                              min(NameVals.QdScaleRange), ...
                                              max(NameVals.QdScaleRange));
fprintf("%s%sSeed: %s\n", INDENT, INDENT, SEED_STR);
fprintf("%s%sNeed OPF success? %s\n", INDENT, INDENT, string(~~NameVals.assertOPFOk));

%% Dispatch building
t_start = tic;
ddata = sdata; if ~isInterIdxMPC(ddata), ddata = ext2int(ddata); end
ddata.Seed = NameVals.Seed;
if ~isempty(NameVals.Seed), rng(NameVals.Seed); end

% Synthesize demand profile
fprintf("%s%sSynthesizing demand profile...", INDENT, INDENT);
N_b = size(ddata.bus, 1);
ddata.bus(:, 3) = ddata.bus(:, 3) .* ( ...
    min(NameVals.PdScaleRange) + (rand([N_b, 1]) * (max(NameVals.PdScaleRange) - min(NameVals.PdScaleRange))));
ddata.bus(:, 4) = ddata.bus(:, 4) .* ( ...
    min(NameVals.QdScaleRange) + (rand([N_b, 1]) * (max(NameVals.QdScaleRange) - min(NameVals.QdScaleRange))));
clear N_b;
fprintf("done.\n");

% Synthesize system parameters
fprintf("%s%sSynthesizing system parameters...", INDENT, INDENT);
N_br = size(ddata.branch, 1);
ddata.branch(:, 3) = ddata.branch(:, 3) .* ( ...
    min(NameVals.RScaleRange) + (rand([N_br, 1]) * (max(NameVals.RScaleRange) - min(NameVals.RScaleRange))));
ddata.branch(:, 4) = ddata.branch(:, 4) .* ( ...
    min(NameVals.XScaleRange) + (rand([N_br, 1]) * (max(NameVals.XScaleRange) - min(NameVals.XScaleRange))));
ddata.branch(:, 5) = ddata.branch(:, 5) .* ( ...
    min(NameVals.BScaleRange) + (rand([N_br, 1]) * (max(NameVals.BScaleRange) - min(NameVals.BScaleRange))));
clear N_br;
fprintf("done.\n");

% Run ACOPF routine
%   just to make sure that the dispatch data is OPF-feasible
fprintf("%s%sRunning ACOPF routine...", INDENT, INDENT);
ddata = runopf(ddata, NameVals.MPOptions);
ddata = ext2int(ddata);
fprintf("done.\n");
if NameVals.assertOPFOk, assert(ddata.success, NO_ACOPF_SOLN_STR); end

%% Summary
fprintf("%sDispatch data built in %f s.\n", INDENT, toc(t_start));

end
