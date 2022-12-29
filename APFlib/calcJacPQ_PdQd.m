function [de_dc] = calcJacPQ_PdQd(Cd)
%CALCJACPQ_PDQD Computes the derivatives of net nodal injections w.r.t. demand draws
%   de_dc = CALCJACPQ_PDQD(Cd)
%   calculates the Jacobian of the vector \( \boldsymbol{e} \) w.r.t. the vector \( \boldsymbol{d} \).
%
%   See Section 2.3.1 of the manuscript.
% 
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7491508
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    Cd {mustBeNonnegative, mustBeInteger}
end

%% Main
[N_b, N_d] = size(Cd);

de_dd = [-Cd, sparse(N_b, N_d); sparse(N_b, N_d), -Cd];

end
