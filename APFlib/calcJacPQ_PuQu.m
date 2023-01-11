function [de_dc] = calcJacPQ_PuQu(Cu)
%CALCJACPQ_PUQU Computes the derivatives of net nodal injections w.r.t. supply injections
%   de_dc = CALCJACPQ_PUQU(Cu)
%   calculates the Jacobian of the vector \( \boldsymbol{e} \) w.r.t. the vector \( \boldsymbol{c} \).
%
%   See Section 2.3.1 of the manuscript.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7523544
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    Cu {mustBeNonnegative, mustBeInteger}
end

%% Main
[N_b, N_u] = size(Cu);

de_dc = [Cu, sparse(N_b, N_u); sparse(N_b, N_u), Cu];

end
