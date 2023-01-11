function [out] = calcMetric(pdata, Vm, Va, Ps, NameVals)
%CALCMETRIC Calculates the dummy metric and its gradient
%   [out] = CALCMETRIC(pdata, Vm, Va, Ps, NameVals)
%
%   Please refer to Section 4.4.1 of the manuscript for more details.
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
    Vm {mustBeNumeric, mustBeVector}
    Va {mustBeNumeric, mustBeVector, mustBeEqualSize(Vm, Va)}
    Ps {mustBeNumeric, mustBeScalarOrEmpty}
    NameVals.getGradient {mustBeNumericOrLogical} = true
    NameVals.ignoreVaRef {mustBeNumericOrLogical} = true
end

%% Main
if isequaln(pdata.Cf, NaN)
    [Cf, Ct] = makeCfCt(pdata.dispatch);
else
    [Cf, Ct] = deal(pdata.Cf, pdata.Ct);
end
Cft = Cf - Ct;
Vm_ = Vm - 1.0;
Va_ = Cft * Va;

out.Fv = norm(Vm_, 2)^2;
out.Fa = norm(Va_, 2)^2;
out.Fs = Ps^2;
out.F = out.Fv + out.Fa + out.Fs;

if NameVals.getGradient
    out.dFv = 2 * Vm_;
    out.dFa = 2 * (Cft.') * Va_;
    out.dFs = 2 * Ps;
    if NameVals.ignoreVaRef
        out.dF = [out.dFv; out.dFa(pdata.are_nonref_buses); out.dFs];
    else
        out.dF = [out.dFv; out.dFa; out.dFs];
    end
else
    [out.dFv, out.dFa, out.dFs, out.dF] = deal(NaN, NaN, NaN, NaN);
end

end
