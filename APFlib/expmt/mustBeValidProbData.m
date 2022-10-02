function mustBeValidProbData(pdata)
%MUSTBEVALIDPROBDATA Validates a set of APF problem data
%   MUSTBEVALIDPROBDATA(pdata)
%   checks the validity of a given problem data struct using `checkProbData`.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~checkProbData(pdata, "beVerbose", false)
    eid = "mustBeValidProbData:notValidProbData";
    msg = "Problem data must be valid.\n";
    msg = sprintf("%sPlease use `checkProbData()` to see the errors.\n", msg);
    throwAsCaller(MException(eid, msg))
end

end
