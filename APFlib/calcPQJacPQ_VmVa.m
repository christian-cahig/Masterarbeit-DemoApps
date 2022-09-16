function [e, J] = calcPQJacPQ_VmVa(vm, va, Y, nref)
%CALCPQJACPQ_VMVA Computes the net nodal injections from and their derivatives w.r.t. bus voltages
%   e = CALCPQJACPQ_VMVA(vm, va, Y)
%   [e, J] = CALCPQJACPQ_VMVA(vm, va, Y)
%   [e, J] = CALCPQJACPQ_VMVA(vm, va, Y, nref)
%   The first output is the vector concatenation of the net active and reactive nodal injections
%   computed using the bus voltage magnitudes, their corresponding phase angles, and the system bus
%   admittance matrix.
%   If the fourth input argument is unspecified, the second output is the corresponding Jacobian of
%   the net active and reactive nodal injections w.r.t. the bus voltage magnitudes and phase angles.
%   If the fourth input argument is specified, the Jacobian is w.r.t. the bus voltage magnitudes
%   and non-reference-bus phase angles.
%
%   See Sections 2.2.2 and 3.2.3 of the manuscript.
% 
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of Anticipatory Power Flow, and is covered by the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    vm {mustBeNumeric, mustBeVector}
    va {mustBeNumeric, mustBeVector, mustBeEqualSize(vm, va)}
    Y {mustBeNumeric}
    nref {mustBePositive, mustBeInteger, mustBeScalarOrEmpty} = []
end

%% Main
vbus = vm .* exp(1j * va);  % Bus voltage phasors
ibus_ = conj(Y * vbus);     % Complex conjugate of net nodal current phasors
sbus = vbus .* ibus_;       % Net complex nodal injections

e = [real(sbus); imag(sbus)];

if nargout > 1
    N_b = numel(vm);
    Vbus = sparse(1:N_b, 1:N_b, vbus, N_b, N_b);
    Ibus_ = sparse(1:N_b, 1:N_b, ibus_, N_b, N_b);
    Vnorm = sparse(1:N_b, 1:N_b, vbus ./ vm, N_b, N_b);

    dsbus_dvm = Vbus * conj(Y * Vnorm) + Ibus_ * Vnorm;
    dsbus_dva = 1j * Vbus * (Ibus_ - conj(Y * Vbus));
    J = [real(dsbus_dvm), real(dsbus_dva); ...
         imag(dsbus_dvm), imag(dsbus_dva)];
    if ~isempty(nref), J(:, N_b + nref) = []; end
end

end
