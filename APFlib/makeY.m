function [Y, Yf, Yt] = makeY(mpc)
%MAKEY Builds the system bus admittance matrix
%   Y = MAKEY(mpc)
%   [Y, Yf, Yt] = MAKEY(mpc)
%   returns the system bus admittance matrix, and, optionally, the from- and to-side system branch
%   admittance matrices, from a given MATPOWER case struct, assumed to be using internal indexing.
%
%   With N_b and N_br as the number of buses and of branches, Y, Yf, and Yt are sparse with dimen-
%   sions of N_b-by-N_b, N_br-by-N_b, and N_br-by-N_b.
% 
%   See Section 2.1.8 of the manuscript for more details.
%
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%      format and indexing, please refer to Appendix B and Section 9.4 of the User's Manual.
%   2. This is a minor modification of MATPOWER's `makeYbus` (Section 9.5.6 of the User's Manual).
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7491508
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeInterIdxMPC}
end

%% Utilities
baseMVA = mpc.baseMVA;
bus     = mpc.bus;
branch  = mpc.branch;
N_b = size(bus, 1);
N_br = size(branch, 1);

% Named indices for the bus and the branch data matrices
[~, ~, ~, ~, ~, ~, ~, ~, GS, BS, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = idx_bus;
[~, ~, BR_R, BR_X, BR_B, ~, ~, ~, ...
    TAP, SHIFT, BR_STATUS, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = idx_brch;

%% Elements of branch admittance matrix
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
Ys = stat ./ (branch(:, BR_R) + 1j * branch(:, BR_X));  %% series admittance
Bc = stat .* branch(:, BR_B);                           %% line charging susceptance
tap = ones(N_br, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
tap = tap .* exp(1j*pi/180 * branch(:, SHIFT)); %% add phase shifters
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

%% Shunt admittances
Ysh = (bus(:, GS) + 1j * bus(:, BS)) / baseMVA; %% vector of shunt admittances

%% Branch connection matrices
[Cf, Ct] = makeCfCt(mpc);

%% System branch admittance matrices
Yf = sparse(1:N_br, 1:N_br, Yff, N_br, N_br) * Cf + sparse(1:N_br, 1:N_br, Yft, N_br, N_br) * Ct;
Yt = sparse(1:N_br, 1:N_br, Ytf, N_br, N_br) * Cf + sparse(1:N_br, 1:N_br, Ytt, N_br, N_br) * Ct;

%% System bus admittance matrix
Y = (Cf' * Yf) +  ...                       % From-side branch admittances
    (Ct' * Yt) +  ...                       % To-side branch admittances
    sparse(1:N_b, 1:N_b, Ysh, N_b, N_b);    % Shunt admittances

end
