function [out] = solveAPFEqns(pdata, confg, NameVals)
%SOLVEAPFEQNS Solves the APF equations
%   out = SOLVEAPFEQNS(pdata, confg, NameVals)
%   solves the APF equations for
%       (i) the voltage magnitudes at all buses,
%       (ii) the voltage phase angles at the non-reference buses, and
%       (iii) the distributed active power slack,
%   using methods provided by the MATLAB Optimization Toolbox.
% 
%   Please refer to Section 3.3 of the manuscript for more details.
% 
%   Notes
%   ---------------------
%   This function requires the MATLAB Optimization Toolbox, and has been developed and tested on
%   version 9.0. In particular, this is a thin wrapper around `fsolve`. For more details, see
%   https://www.mathworks.com/help/releases/R2021b/optim/ug/fsolve.html.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7180586
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    pdata
    confg
    NameVals.Pu {mustBeNumeric, mustBeVector} = pdata.snapshot.gen(:, 2) ./ pdata.sys_baseMVA
    NameVals.Qu {mustBeNumeric, mustBeVector} = pdata.snapshot.gen(:, 2) ./ pdata.sys_baseMVA
    NameVals.PsDistribs {mustBeNonnegative, mustBeVector} = ( ...
        (pdata.snapshot.gen(:, 9) - pdata.snapshot.gen(:, 10)) ...
        / sum(pdata.snapshot.gen(:, 9) - pdata.snapshot.gen(:, 10)) ) / pdata.sys_baseMVA
    NameVals.InitVm {mustBePositive, mustBeVector} = pdata.snapshot.bus(:, 8)
    NameVals.InitVa {mustBeNumeric, mustBeVector} = deg2rad(pdata.snapshot.bus(:, 9))
    NameVals.InitPs (1, 1) {mustBeNumeric} = 0
end

% Check vector sizes
if numel(NameVals.Pu) ~= pdata.N_u
    msg = 'Pu'; N = pdata.N_u;
elseif numel(NameVals.Qu) ~= pdata.N_u
    msg = 'Qu'; N = pdata.N_u;
elseif numel(NameVals.PsDistribs) ~= pdata.N_u
    msg = 'PsDistribs'; N = pdata.N_u;
elseif numel(NameVals.InitVm) ~= pdata.N_b
    msg = 'InitVm'; N = pdata.N_b;
elseif numel(NameVals.InitVa) ~= pdata.N_b
    msg = 'InitVa'; N = pdata.N_b;
else
    msg = '';
end
if ~isempty(msg)
    msg = sprintf("Invalid name-value argument %s.", msg);
    msg = sprintf("%s Value must be an %d-vector.\n", msg, N);
    throw(MException("solveAPFEqns:inputError", msg));
end
clear msg;

% Check validity
if abs(sum(NameVals.PsDistribs)) - 1.0 > 1e-7
    msg = "Invalid name-value argument 'PsDistribs'.";
    msg = sprintf("%s abs(sum(PsDistribs) - 1.0) > 1e-7.\n", msg);
    throw(MException("solveAPFEqns:inputError", msg));
end
out.PsDistribs = NameVals.PsDistribs;

%% Intro
% Net nodal injections according to supply injections and demand draws
if isequaln(pdata.Cu, NaN), Cu = makeCu(pdata.dispatch); else, Cu = pdata.Cu; end
if isequaln(pdata.Cd, NaN), Cd = makeCd(pdata.dispatch); else, Cd = pdata.Cd; end
e = calcPQFromCD(Cu, Cd, NameVals.Pu, NameVals.Qu, ...
                 pdata.dispatch.bus(pdata.dem_buses, 3) / pdata.sys_baseMVA, ...
                 pdata.dispatch.bus(pdata.dem_buses, 4) / pdata.sys_baseMVA);
clear Cd;

% Per-bus slack distribution
CuKappa = (-Cu) * out.PsDistribs;
clear Cu;

% System bus admittance matrix and non-reference buses
if isequaln(pdata.Y, NaN), Y = makeY(pdata.dispatch); else, Y = pdata.Y; end
are_nonRef_buses = (1 : pdata.N_b) ~= pdata.ref_bus;

% Initial iterates
[out.Vm, out.Va, out.Ps] = deal(NameVals.InitVm, NameVals.InitVa, NameVals.InitPs);
out.Va = out.Va + (pdata.va_ref - out.Va(pdata.ref_bus));

%% Solve
Phi = @(x)evalAPFEqns(x, e, Y, pdata.ref_bus, pdata.va_ref, CuKappa, pdata.N_b);
t = tic;
[x, out.Phi, out.Status, out_, out.Jac] = fsolve( ...
    Phi, [out.Vm; out.Va(are_nonRef_buses); out.Ps], confg);
out.Time = toc(t);
out.Vm = x(1 : pdata.N_b);
out.Va(are_nonRef_buses) = x(pdata.N_b + 1 : (2 * pdata.N_b) - 1);
out.Ps = x(end);

%% Summary
if norm(out.Phi, Inf) <= sqrt(confg.FunctionTolerance)
    out.converged = true;
elseif (out.Status < 5) && (norm(out.Phi, 2)^2 <= confg.FunctionTolerance)
    out.converged = true;
else
    out.converged = false;
end
out.Steps = out_.iterations;
out.FuncEvals = out_.funcCount;
out.Message = out_.message;

end
