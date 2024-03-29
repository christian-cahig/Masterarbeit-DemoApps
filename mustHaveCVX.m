function mustHaveMATPOWER(~)
%MUSTHAVECVX Validates CVX installation
%   MUSTHAVECVX()
%   checks if CVX is installed via looking for `cvx_version` in the MATLAB search path.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7523544
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~exist('cvx_version', 'file')
    eid = "mustHaveCVX:noCVX";
    msg = "CVX must be installed first.\n";
    throwAsCaller(MException(eid, msg));
end

end
