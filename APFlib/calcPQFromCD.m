function [p, q] = calcPQFromCD(Cu, Cd, pu, qu, pd, qd)
%CALCEFROMCD Calculates the net nodal injections from supply injections and demand draws
%   [p, q] = CALCEFROMCD(Cu, Cd, pu, qu, pd, qd)
%   returns the net active and reactive nodal injections.
%
%   [e] = CALCEFROMCD(Cu, Cd, pu, qu, pd, qd)
%   returns the concatenation of the vectors of net active and reactive nodal injections.
%
%   Please refer to Section 2.2.1 of the manuscript.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    Cu {mustBeNonnegative, mustBeInteger}
    Cd {mustBeNonnegative, mustBeInteger}
    pu {mustBeNumeric, mustBeVector}
    qu {mustBeNumeric, mustBeVector, mustBeEqualSize(pu, qu)}
    pd {mustBeNumeric, mustBeVector}
    qd {mustBeNumeric, mustBeVector, mustBeEqualSize(pd, qd)}
end

%% Main
if nargout == 1
    p = [(Cu * pu) - (Cd * pd); ...
         (Cu * qu) - (Cd * qd)];
else
    p = (Cu * pu) - (Cd * pd);
    q = (Cu * qu) - (Cd * qd);
end

end
