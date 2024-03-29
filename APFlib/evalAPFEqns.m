function [phi, J] = evalAPFEqns(x, e, Y, ref_bus, va_ref, ps_ratio, N_b)
%EVALAPFEQNS Evaluates the APF equations and optionally their Jacobians
% 
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7523544
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Bus voltage phase angles
va = [x(N_b + 1 : N_b + ref_bus - 1);
      va_ref;
      x(N_b + ref_bus : (2 * N_b) - 1)];

%% Nodal power residuals and, optionally, their Jacobians
if nargout > 1
    % (↓) Approach #1
    [phi, J] = calcPQJacPQ_VmVa(x(1 : N_b), va, Y, ref_bus);
    J = [J, [ps_ratio; sparse(N_b, 1)]];
    % (↓) Approach #2
    % [phi, J] = calcPQJacPQ_VmVa(x(1 : N_b), va, Y);
    % J = [J(:, 1 : N_b + ref_bus - 1), ...
    %      J(:, N_b + ref_bus + 1 : end), ...
    %      [ps_ratio; sparse(N_b, 1)]];
else
    phi = calcPQJacPQ_VmVa(x(1 : N_b), va, Y);
end
phi(1 : N_b) = phi(1 : N_b) + (x(end) * ps_ratio);
phi = phi - e;

end
