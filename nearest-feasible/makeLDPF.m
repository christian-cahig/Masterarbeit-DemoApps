function [C1, C2, C3, C4] = makeLDPF(Y, ref_bus)
%MAKELDPF Calculates the linearized decoupled power flow matrices
%   [C1, C2, C3, C4] = MAKELDPF(Y)
%   calculates the linearized decoupled power flow (LDPF) matrices from the bus admittance matrix.
%   If pnet + j*qnet is the vector of net complex power injections at the N_b buses, and the bus
%   voltage phasors have magnitudes vm and phase angle va, then the LDPF model is
%       [pnet; qnet] = [C1, C2; C3, C4] [va; vm]
%   implying that the four submatrices have dimensions N_b-by-N_b.
% 
%   [C1r, C2r, C3r, C4] = MAKELDPF(Y, rbus)
%   returns the reduced LDPF matrices which are formed as follows:
%       C1r = C1; C1r(rbus, :) = []; C1r(:, rbus) = [];
%       C2r = C2; C2r(rbus, :) = [];
%       C3r = C3; C3r(:, rbus) = [];
%   This corresponds to the reduced LDPF model where the net active power injection and the voltage
%   phase angle at the reference bus rbus are excluded.
% 
%   [C] = MAKELDPF(...)
%   returns either [C1, C2; C3, C4] or [C1r, C2r; C3r, C4].
% 
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    Y {mustBeNumeric}
    ref_bus {mustBePositive, mustBeInteger, mustBeScalarOrEmpty} = []
end

%% Main
N_b = size(Y, 1);
[G, B] = deal(real(Y), imag(Y));
G0 = spdiags(sparse(N_b, 1), 0, G);
B0 = spdiags(sparse(N_b, 1), 0, B);
clear N_b;

[C1, C2, C4] = deal(-B0, -G0, B0); C3 = C2;
C1 = spdiags(sum(B0, 2), 0, C1); clear B0;
C3 = spdiags(sum(G0, 2), 0, C3); clear G0;
C2 = spdiags(sum(G, 2), 0, C2); clear G;
C4 = spdiags(-sum(B, 2), 0, C4); clear B;

if ~isempty(ref_bus)
    C1(ref_bus, :) = []; C1(:, ref_bus) = [];
    C2(ref_bus, :) = []; C3(:, ref_bus) = [];
end

if nargout < 4, C1 = [C1, C2; C3, C4]; end

end
