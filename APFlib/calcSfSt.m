function [sf, st] = calcSfSt(vm, va, Yf, Yt, fbus, tbus)
%CALCSFST Calculates the branch complex power flows
%   [sf, st] = CALCSFST(vm, va, Yf, Yt, fbus, tbus)
%   returns the vectors of the from- and of the to-side complex branch power flows computed using
%   the bus voltage magnitudes, their corresponding phase angles and the from- and the to-side
%   system branch admittance matrices.
% 
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format and indexing, please refer to Appendix B and Section 9.4.
%   2. This is a tweaked version of MATPOWER's `opf_branch_flow_fcn` function.
% 
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7180586
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    vm {mustBeNumeric, mustBeVector}
    va {mustBeNumeric, mustBeVector, mustBeEqualSize(vm, va)}
    Yf {mustBeNumeric}
    Yt {mustBeNumeric, mustBeEqualSize(Yf, Yt)}
    fbus {mustBeInteger, mustBeVector}
    tbus {mustBeInteger, mustBeVector, mustBeEqualSize(fbus, tbus)}
end

%% Main
v = vm .* exp(1j * va);
sf = v(fbus) .* conj(Yf * v);
st = v(tbus) .* conj(Yt * v);

end
