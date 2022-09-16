function [out] = solveExtdEconDisp(pdata, NameVals)
%SOLVEEXTDECONDISP Anticipates supply injections via the extended economic dispatch
%   [out] = SOLVEEXTDECONDISP(pdata, NameVals)
%   anticipates the supply injections by solving an extended variant of economic dispatch.
% 
%   Please refer to Section 3.2 of the manuscript for more details.
% 
%   Notes
%   ---------------------
%   This function relies on CVX, and has been developed and tested on version 2.2, build 1148.
%   The accompanying User's Guide is freely available at
%   http://web.cvxr.com/cvx/doc/CVX.pdf.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
% 
%   This file is part of Anticipatory Power Flow, and is covered by the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    pdata
    NameVals.Ph {mustBeNonnegative, mustBeScalarOrEmpty} = []
    NameVals.Qh {mustBeNumeric, mustBeScalarOrEmpty} = []
    NameVals.Po {mustBePositive, mustBeScalarOrEmpty} = []
    NameVals.Qo {mustBeNumeric, mustBeScalarOrEmpty} = []
    NameVals.PuReg {mustBeNonnegative} = 1e-3
    NameVals.QuReg {mustBeNonnegative} = 1e-3
    NameVals.useCQP (1, 1) {mustBeNumericOrLogical} = false
end

%% Intro
[out.PuReg, out.QuReg] = deal(NameVals.PuReg, NameVals.QuReg);

% Anticipated shunt consumption
if isempty(NameVals.Ph) && isempty(NameVals.Qh)
    [out.Ph, out.Qh] = calcPhQh(pdata.dispatch, "Vm", pdata.snapshot.bus(:, 8));
elseif isempty(NameVals.Ph)
    [out.Ph, ~] = calcPhQh(pdata.dispatch, "Vm", pdata.snapshot.bus(:, 8));
    out.Qh = NameVals.Qh;
elseif isempty(NameVals.Qh)
    out.Ph = NameVals.Ph;
    [~, out.Qh] = calcPhQh(pdata.dispatch, "Vm", pdata.snapshot.bus(:, 8));
else
    [out.Ph, out.Qh] = deal(NameVals.Ph, NameVals.Qh);
end

% Anticipated system losses
if isempty(NameVals.Po) && isempty(NameVals.Qo)
    [out.Po, out.Qo] = calcLoss(pdata.dispatch, ...
        "Vm", pdata.snapshot.bus(:, 8), "Va", pdata.snapshot.bus(:, 9));
elseif isempty(NameVals.Po)
    [out.Po, ~] = calcLoss(pdata.dispatch, ...
        "Vm", pdata.snapshot.bus(:, 8), "Va", pdata.snapshot.bus(:, 9));
    out.Qo = NameVals.Qo;
elseif isempty(NameVals.Qo)
    out.Po = NameVals.Po;
    [~, out.Qo] = calcLoss(pdata.dispatch, ...
        "Vm", pdata.snapshot.bus(:, 8), "Va", pdata.snapshot.bus(:, 9));
else
    [out.Po, out.Qo] = deal(NameVals.Po, NameVals.Qo);
end

% Parameters
Pu_snap = pdata.snapshot.gen(:, 2) ./ pdata.sys_baseMVA;
Qu_snap = pdata.snapshot.gen(:, 3) ./ pdata.sys_baseMVA;

Pu_min = pdata.dispatch.gen(:, 10) ./ pdata.sys_baseMVA;
nlbp_units = find(Pu_min == -Inf);
lbp_units = find(Pu_min ~= -Inf);
Pu_min(nlbp_units) = min(Pu_snap(nlbp_units), min(Pu_min(lbp_units)));
out.Pu_min = Pu_min;
clear nlbp_units lbp_units Pu_min;

Pu_max = pdata.dispatch.gen(:, 9) ./ pdata.sys_baseMVA;
nubp_units = find(Pu_max == Inf);
ubp_units = find(Pu_max ~= Inf);
Pu_max(nubp_units) = max(Pu_snap(nubp_units), max(Pu_max(ubp_units)));
out.Pu_max = Pu_max;
clear nubp_units ubp_units Pu_max;

Qu_min = pdata.dispatch.gen(:, 5) ./ pdata.sys_baseMVA;
nlbq_units = find(Qu_min == -Inf);
lbq_units = find(Qu_min ~= -Inf);
Qu_min(nlbq_units) = min(Qu_snap(nlbq_units), min(Qu_min(lbq_units)));
out.Qu_min = Qu_min;
clear nlbq_units lbq_units Qu_min;

Qu_max = pdata.dispatch.gen(:, 4) ./ pdata.sys_baseMVA;
nubq_units = find(Qu_max == Inf);
ubq_units = find(Qu_max ~= Inf);
Qu_max(nubq_units) = max(Qu_snap(nubq_units), max(Qu_max(ubq_units)));
out.Qu_max = Qu_max;
clear nubq_units ubq_units Qu_max;

out.p_need = sum(pdata.dispatch.bus(pdata.dem_buses, 3) ./ pdata.sys_baseMVA) ...
    + out.Ph + out.Po;
out.p_need = max(sum(out.Pu_min), min(sum(out.Pu_max), out.p_need));
out.q_need = sum(pdata.dispatch.bus(pdata.dem_buses, 4) ./ pdata.sys_baseMVA) ...
    + out.Qh + out.Qo;
out.q_need = max(sum(out.Qu_min), min(sum(out.Qu_max), out.q_need));

%% Solve
if NameVals.useCQP
    A = [ones(1, pdata.N_u),    sparse(1, pdata.N_u); ...
         sparse(1, pdata.N_u),  ones(1, pdata.N_u)];
    a = [p_need; q_need];
    mu = zeros([2*pdata.N_u, 1]);
    mu(1 : pdata.N_u) = sqrt(out.PuReg);
    mu(pdata.N_u + 1 : end) = sqrt(out.QuReg);
    H = sparse([1 : 2*pdata.N_u], [1 : 2*pdata.N_u], 1 + mu, 2*pdata.N_u, 2*pdata.N_u);
    h = -2 * (mu .* [Pu_snap; Qu_snap]);
    clear mu;

    t0 = tic;
    cvx_begin quiet
        variable w(2*pdata.N_u, 1)
        minimize( quad_form(w, H) + sum(h .* w) );
        subject to
            A * w == a;
            [out.Pu_min; out.Qu_min] <= w <= [out.Pu_max; out.Qu_max];
    cvx_end
else
    t0 = tic;
    cvx_begin quiet
        variable Pu(pdata.N_u, 1)
        variable Qu(pdata.N_u, 1)
        minimize( norm(Pu, 2) + norm(Qu, 2) ...
            + (out.PuReg * norm(Pu - Pu_snap, 2)) ...
            + (out.QuReg * norm(Qu - Qu_snap, 2)) );
        subject to
            sum(Pu) == out.p_need;
            sum(Qu) == out.q_need;
            out.Pu_min <= Pu <= out.Pu_max;
            out.Qu_min <= Qu <= out.Qu_max;
    cvx_end
end

%% Summary
out.Time = toc(t0);
if NameVals.useCQP
    [out.Pu, out.Qu] = deal(w(1 : pdata.N_u), w(pdata.N_u + 1 : end));
else
    [out.Pu, out.Qu] = deal(Pu, Qu);
end
out.Status = cvx_status;
out.Steps = cvx_slvitr;
out.OptVal = cvx_optval;

end
