function mustBeMPC(mpc)
%MUSTBEMPC Validates a MATPOWER case
%   MUSTBEMPC(mpc) checks if `mpc` is a valid MATPOWER case struct using `isMPC`.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7523544
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~isMPC(mpc)
    eid = "mustBeMPC:notMPC";
    msg = "Must be a MATPOWER case struct.\n";
    throwAsCaller(MException(eid, msg))
end

end
