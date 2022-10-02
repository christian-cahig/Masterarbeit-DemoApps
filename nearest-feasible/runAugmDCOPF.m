function [out] = runAugmDCOPF(pdata, NameVals)
%RUNAUGMDCOPF Fancy wrapper for `rundcopf`
%   [out] = RUNAUGMDCOPF(pdata, NameVals)
% 
%   Please refer to Section 4.3.2 of the manuscript for more details.
% 
%   Notes
%   ---------------------
%   This function requires MATPOWER. It has been developed and tested on version 7.1.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
% 
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    pdata
    NameVals.MPOptions struct = getGoToMPOpts()
    NameVals.QuStrat (1, 1) {mustBeInteger} = 1
end

%% Main
mpc = rundcopf(pdata.dispatch, NameVals.MPOptions);
out.Pu = mpc.var.val.Pg;
out.Va = mpc.var.val.Va;
NameVals.MPOptions.pf.alg = 'FDBX';
NameVals.MPOptions.pf.fd.max_it = 100;
mpc_ = runpf(mpc, NameVals.MPOptions);
if mpc_.success
    out.Vm = mpc_.bus(:, 8);
    % out.Va = deg2rad(mpc_.bus(:, 9));
    % out.Pu = mpc_.gen(:, 2) ./ pdata.sys_baseMVA;
    out.Qu = mpc_.gen(:, 3) ./ pdata.sys_baseMVA;
else
    % out.Vm = ones([pdata.N_b, 1]);
    % out.Va = mpc.var.val.Va;
    % out.Pu = mpc.var.val.Pg;
    if NameVals.QuStrat == 0
        out.Qu = pdata.snapshot.gen(:, 3) ./ pdata.sys_baseMVA;
    elseif NameVals.QuStrat == 1
        tan_pf_ang = sum(pdata.snapshot.gen(:, 3)) ./ sum(pdata.snapshot.gen(:, 2));
        out.Qu = out.Pu .* tan_pf_ang;
    elseif NameVals.QuStrat == 2
        tan_pf_ang = sum(pdata.dispatch.bus(:, 4)) ./ sum(pdata.dispatch.bus(:, 3));
        out.Qu = out.Pu .* tan_pf_ang;
    else
        out.Qu = pdata.snapshot.gen(:, 3) ./ pdata.sys_baseMVA;
    end
end
out.SPFOk = mpc_.success;

end
