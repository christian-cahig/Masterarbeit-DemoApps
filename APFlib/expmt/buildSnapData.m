function [sdata] = buildSnapData(mpc, NameVals)
%BUILDSNAPDATA Derives snapshot data from a MATPOWER case
%   sdata = BUILDSNAPDATA(mpc, NameVals)
%   synthesizes a set of snapshot data based on a given MATPOWER case.
%
%   Inputs
%   ---------------------
%   mpc : the MATPOWER case struct forming the basis for the snapshot data.
%   NameVals : optional name-value pairs.
%   (-) MPOptions, a struct containing the MATPOWER options. If unspecified, the one provided by
%       `getGoToMPOpts` is used.
%   (-) RScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.90; 1.10]` by default) for the branch series resistances in `mpc`.
%   (-) XScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.90; 1.10]` by default) for the branch series reactances in `mpc`.
%   (-) BScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.95; 1.05]` by default) for the branch shunt reactances in `mpc`.
%   (-) PdScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.80; 1.20]` by default) for the active demand draws in `mpc`.
%   (-) PdShiftRange, a 2-element positive vector containing the minimum and the maximum offsets
%       (`[-1e-4; 1e-4]` by default) for the per-unit active demand draws in `mpc`.
%   (-) QdScaleRange, a 2-element positive vector containing the minimum and the maximum multipliers
%       (`[0.80; 1.20]` by default) for the reactive demand draws in `mpc`.
%   (-) QdShiftRange, a 2-element positive vector containing the minimum and the maximum offsets
%       (`[-1e-4; 1e-4]` by default) for the per-unit reactive demand draws in `mpc`.
%   (-) Seed, a nonnegative integer with which the random number generator is seeded. If unspecified
%       the random number generator is not seeded.
%
%   Outputs
%   ---------------------
%   sdata : a MATPOWER case struct (in MATPOWER's internal indexing) containing the snapshot data.
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
%   https://doi.org/10.5281/zenodo.7491508
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeMPC} = case14()
    NameVals.MPOptions struct = getGoToMPOpts()
    NameVals.RScaleRange (2, 1) {mustBePositive} = [0.90; 1.10]
    NameVals.XScaleRange (2, 1) {mustBePositive} = [0.90; 1.10]
    NameVals.BScaleRange (2, 1) {mustBePositive} = [0.95; 1.05]
    NameVals.PdScaleRange (2, 1) {mustBePositive} = [0.80; 1.20]
    NameVals.PdShiftRange (2, 1) {mustBeNumeric} = [-1e-4; 1e-4]
    NameVals.QdScaleRange (2, 1) {mustBePositive} = [0.80; 1.20]
    NameVals.QdShiftRange (2, 1) {mustBeNumeric} = [-1e-4; 1e-4]
    NameVals.Seed {mustBeInteger, mustBeNonnegative, mustBeScalarOrEmpty} = []
end

%% Utilities
INDENT = "   ";
TRY_AGAIN_STR = "Try running again (if unseeded) or use a different seed.";
NO_ACOPF_SOLN_STR = sprintf("%s%sNo ACOPF solution found.\n%s%s\n", INDENT, INDENT, INDENT, TRY_AGAIN_STR);
if isempty(NameVals.Seed), SEED_STR = "None"; else, SEED_STR = sprintf("%d", NameVals.Seed); end

%% Intro
fprintf("%sBuilding snapshot data from a MATPOWER case.\n", INDENT);
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
fprintf("%s%sDemand P shift: (%.3e, %.3e)\n", INDENT, INDENT, ...
                                              min(NameVals.PdShiftRange), ...
                                              max(NameVals.PdShiftRange));
fprintf("%s%sDemand Q scale: (%.3f, %.3f)\n", INDENT, INDENT, ...
                                              min(NameVals.QdScaleRange), ...
                                              max(NameVals.QdScaleRange));
fprintf("%s%sDemand Q shift: (%.3e, %.3e)\n", INDENT, INDENT, ...
                                              min(NameVals.QdShiftRange), ...
                                              max(NameVals.QdShiftRange));
fprintf("%s%sSeed: %s\n", INDENT, INDENT, SEED_STR);

%% Snapshot building
t_start = tic;
sdata = ext2int(mpc); sdata.Seed = NameVals.Seed;
if ~isempty(NameVals.Seed), rng(NameVals.Seed); end

% Prepare branch data
fprintf("%s%sPreparing branch data...", INDENT, INDENT);
N_br = size(sdata.branch, 1);

nu = min(NameVals.RScaleRange) + rand(N_br, 1) * (          ...
    max(NameVals.RScaleRange) - min(NameVals.RScaleRange)   ...
);
sdata.branch(:, 3) = nu .* sdata.branch(:, 3);

iota = min(NameVals.XScaleRange) + rand(N_br, 1) * (        ...
    max(NameVals.XScaleRange) - min(NameVals.XScaleRange)   ...
);
sdata.branch(:, 4) = iota .* sdata.branch(:, 4);

tau = min(NameVals.BScaleRange) + rand(N_br, 1) * (         ...
    max(NameVals.BScaleRange) - min(NameVals.BScaleRange)   ...
);
sdata.branch(:, 5) = tau .* sdata.branch(:, 5);

clear N_br nu iota tau;
fprintf("done.\n");

% Prepare demand profile
fprintf("%s%sPreparing demand profile...", INDENT, INDENT);
[dem_buses, N_d] = getDemBuses(sdata);

alpha_mul = min(NameVals.PdScaleRange) + rand(N_d, 1) * (   ...
    max(NameVals.PdScaleRange) - min(NameVals.PdScaleRange) ...
);
sdata.bus(dem_buses, 3) = alpha_mul .* sdata.bus(dem_buses, 3);
alpha_off = min(NameVals.PdShiftRange) + rand(N_d, 1) * (   ...
    max(NameVals.PdShiftRange) - min(NameVals.PdShiftRange) ...
);
sdata.bus(dem_buses, 3) = sdata.bus(dem_buses, 3) + (sdata.baseMVA * alpha_off);

omega_mul = min(NameVals.QdScaleRange) + rand(N_d, 1) * (   ...
    max(NameVals.QdScaleRange) - min(NameVals.QdScaleRange)  ...
);
sdata.bus(dem_buses, 4) = omega_mul .* sdata.bus(dem_buses, 4);
omega_off = min(NameVals.QdShiftRange) + rand(N_d, 1) * (   ...
    max(NameVals.QdShiftRange) - min(NameVals.QdShiftRange) ...
);
sdata.bus(dem_buses, 3) = sdata.bus(dem_buses, 3) + (sdata.baseMVA * omega_off);

clear dem_buses N_d alpha_mul alpha_off omega_mul omega_off;
fprintf("done.\n");

% Run ACOPF routine
fprintf("%s%sRunning ACOPF routine...", INDENT, INDENT);
sdata = runopf(sdata, NameVals.MPOptions);
sdata = ext2int(sdata);
fprintf("done.\n");
assert(sdata.success, NO_ACOPF_SOLN_STR);

%% Summary
fprintf("%sSnapshot data built in %f s.\n", INDENT, toc(t_start));

end
