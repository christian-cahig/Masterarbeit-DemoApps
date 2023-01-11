function [po, qo] = calcPoQo(mpc, NameVals)
%CALCPOQO Calculates the total power lost in the branches
%   [po, qo] = CALCPOQO(mpc)
%   returns the total (per-unit) active and reactive losses in the branches of a MATPOWER case.
% 
%   [po, qo] = CALCPOQO(mpc, "Vm", vm, "Va", va)
%   calculates the total losses from a user-defined set of bus voltage magnitudes, of bus voltage
%   phase angles, or both, but using the branch parameters in a given (compatible) MATPOWER case.
%   This is equivalent to the logic
%       mpc.bus(:, 8) = vm;
%       mpc.bus(:, 9) = va;
%       [po, qo] = calcPoQo(mpc);
%   which implies that the given voltage phase angles are assumed to be in degrees. This syntax is
%   intended for evaluating losses at different bus voltage profiles but fixed system parameters.
% 
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%      format, please refer to Appendix B of the User's Manual.
%   2. This is a thin wrapper around MATPOWER's built-in utility `get_losses` (see Section 9.2.4 of
%      the User's Manual).
% 
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7523544
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeMPC}
    NameVals.Vm {mustBePositive, mustBeVector} = []
    NameVals.Va {mustBeNumeric, mustBeVector} = []
end

if ~isInterIdxMPC(mpc), mpc = ext2int(mpc); end;
N_b = size(mpc.bus, 1);
if ~isempty(NameVals.Vm) && (numel(NameVals.Vm) ~= N_b)
    NON_NB_VECTOR = "Vm";
elseif ~isempty(NameVals.Va) && (numel(NameVals.Va) ~= N_b)
    NON_NB_VECTOR = "Va";
else
    NON_NB_VECTOR = '';
end
if ~isempty(NON_NB_VECTOR)
    msg = sprintf("Invalid name-value argument %s.", NON_NB_VECTOR);
    msg = sprintf("%s Value must be a %d-vector.\n", msg, N_b);
    throw(MException("calcPoQo:inputError", msg));
end
clear N_b NON_NB_VECTOR;

%% Main
bus = mpc.bus;
if ~isempty(NameVals.Vm), bus(:, 8) = NameVals.Vm; end;
if ~isempty(NameVals.Va), bus(:, 9) = NameVals.Va; end
[s, c] = get_losses(mpc.baseMVA, bus, mpc.branch);
clear bus;
s = s / mpc.baseMVA;
c = c / mpc.baseMVA;

po = sum(real(s));
qo = sum(imag(s)) - sum(c);

end
