function [out] = makeGSDF(ddata, NameVals)
%MAKEGSDF Calculates the generation shift distribution factors
%   [out] = MAKEGSDF(ddata)
%   calculates the generation shift distribution factors (GSDFs) of a dispatch data set given as a
%   MATPOWER case. The output is a struct with four fields:
%   (-) pp, a matrix collecting the GSDFs from the total active supply injections at the buses to
%       the active power flow across the branches.
%   (-) pq, a matrix collecting the GSDFs from the total reactive supply injections at the buses to
%       the active power flow across the branches.
%   (-) qp, a matrix collecting the GSDFs from the total active supply injections at the buses to
%       the reactive power flow across the branches.
%   (-) qq, a matrix collecting the GSDFs from the total reactive supply injections at the buses to
%       the reactive power flow across the branches.
%   With N_b and N_br as the number of buses and of branches, the GSDF matrices are N_br-by-N_b.
% 
%   [out] = MAKEGSDF(ddata, "Ybus", Y, "RefBus", rbus)
%   skips the calculation of the system bus admittance matrix from scratch, and, optionally, use a
%   reference bus different from that indicated in the MATPOWER case.
% 
%   [out] = MAKEGSDF(..., "ignoreQft", true)
%   skips the calculation of the GSDFs from total supply injections at the buses to the reactive
%   power flows across the branches. Consequently, out.qp and out.qq are set to NaN.
% 
%   The GSDF calculation here is a slight modification of that used in Section II-C of:
%       Q. Hou, N. Zhang, J. Yang, C. Kang, Q. Xia, and M. Miao
%       "Linearized Model for Active and Reactive LMP Considering Bus Voltage Constraints"
%       Proceedings of the 2018 IEEE Power & Energy Society General Meeting (PESGM)
%       DOI: 10.1109/pesgm.2018.8586191
% 
%   Notes
%   ---------------------
%   This function relies heavily on the MATPOWER package. It has been developed and tested on
%   version 7.1. The accompanying User's Manual is freely available at
%   https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format, please refer to Appendix B of the User's Manual.
% 
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    ddata
    NameVals.Y {mustBeNumeric} = NaN
    NameVals.RefBus {mustBePositive, mustBeInteger} = 1
    NameVals.ignoreQft (1, 1) {mustBeNumericOrLogical} = true
    NameVals.getX (1, 1) {mustBeNumericOrLogical} = false
end

if ~isInterIdxMPC(ddata), ddata = ext2int(ddata); end;
if isequaln(NameVals.Y, NaN), NameVals.Y = makeY(ddata); end;

%% Main
% LPDF matrix
C = makeLDPF(NameVals.Y, NameVals.RefBus);

% This is Equation (10) of the Hou et al. paper
N_b = size(ddata.bus, 1);
out.X = sparse(2*N_b, 2*N_b);
if isempty(NameVals.RefBus)
    idxs = (1 : 2*N_b); 
else
    idxs = [find((1 : N_b) ~= NameVals.RefBus), (N_b + 1 : 2*N_b)];
end
out.X(idxs, idxs) = -mldivide(C, speye(size(C, 1)));
if strcmp(lastwarn, 'Matrix is singular to working precision.')
    [out.X, out.pp, out.pq, out.qp, out.qq] = deal(NaN);
    out.success = false;
    return;
end
out.success = true;
clear C;

% GSDF auxiliaries
N_br = size(ddata.branch, 1); branches = (1 : N_br);
z = (ddata.branch(:, 3).^2) + (ddata.branch(:, 4).^2);
diag_g = sparse(branches, branches, ddata.branch(:, 3) ./ z, N_br, N_br);
diag_b = sparse(branches, branches, ddata.branch(:, 4) ./ z, N_br, N_br);
clear z;
[Cf, Ct] = makeCfCt(ddata); Cft = Cf - Ct;
clear Cf Ct;
G1 = Cft * out.X(1:N_b, 1:N_b);
G2 = Cft * out.X(1:N_b, N_b + 1 : end);
G3 = Cft * out.X(N_b + 1 : end, 1:N_b);
G4 = Cft * out.X(N_b + 1 : end, N_b + 1 : end);

% GSDF p-to-p, see Equation (16) of the Hou et al. paper
out.pp = (diag_g * G3) + (diag_b * G1);

% GSDF q-to-p, see Equation (17) of the Hou et al. paper
out.pq = (diag_g * G4) + (diag_b * G2);

if ~NameVals.getX, out.X = NaN; end

if ~NameVals.ignoreQft
    % GSDF p-to-q, see Equation (18) of the Hou et al. paper
    out.qp = (-diag_g * G1) + (diag_b * G3);

    % GSDF q-to-q, see Equation (19) of the Hou et al. paper
    out.qq = (-diag_g * G2) + (diag_b * G4); 
else
    [out.qp, out.qq] = deal(NaN);
end

end
