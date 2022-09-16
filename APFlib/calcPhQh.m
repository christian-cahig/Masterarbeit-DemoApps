function [ph, qh] = calcPhQh(mpc, NameVals)
%CALCPHQH Calculates the total power consumed by the shunt elements
%   [ph, qh] = CALCPHQH(mpc)
%   returns the total (per-unit) active and reactive powers consumed by the shunt elements of a
%   MATPOWER case.
% 
%   [ph, qh] = CALCPHQH(mpc, "Vm", vm, "Va", va)
%   calculates the total shunt-element consumption from a given set of bus voltage magnitudes but
%   using the shunt admittances in a given (compatible) MATPOWER case. This is equivalent to
%       mpc.bus(:, 8) = vm;
%       [ph, qh] = calcPhQh(mpc);
%   This syntax is intended for evaluating shunt consumptions at different bus voltage profiles but
%   fixed shunt-element parameters.
% 
%   Notes
%   ---------------------
%   This function relies heavily on the MATPOWER package. It has been developed and tested on
%   version 7.1. The accompanying User's Manual is freely available at
%   https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format, please refer to Appendix B of the User's Manual. For more details on how MATPOWER models
%   shunt elements, please refer to Section 3.5 of the User's Manual.
% 
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of Anticipatory Power Flow, and is covered by the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeMPC}
    NameVals.Vm {mustBePositive, mustBeVector} = []
end

if ~isInterIdxMPC(mpc), mpc = ext2int(mpc); end;
if ~isempty(NameVals.Vm) && (numel(NameVals.Vm) ~= size(mpc.bus, 1))
    msg = sprintf("Invalid name-value argument Vm.");
    msg = sprintf("%s Value must be a positive %d-vector.\n", msg, size(mpc.bus, 1));
    throw(MException("calcPhQh:inputError", msg));
end

%% Main
bus = mpc.bus;
if ~isempty(NameVals.Vm), bus(:, 8) = NameVals.Vm; end;

Vm_ = mpc.bus(:, 8).^2;
ph = sum(Vm_ .* (mpc.bus(:, 5) / mpc.baseMVA), 1);
qh = -sum(Vm_ .* (mpc.bus(:, 6) / mpc.baseMVA), 1);

end
