function mustBeValidDispData(ddata)
%MUSTBEVALIDDISPDATA Validates a set of dispatch data
%   MUSTBEVALIDDISPDATA(ddata)
%   checks the validity of a given dispatch data struct using `checkDispData`.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~checkDispData(ddata, "beVerbose", false)
    eid = "mustBeValidDispData:notValidDispData";
    msg = "Dispatch data must be valid.\n";
    msg = sprintf("%sPlease use `checkDispData` to see the errors.\n", msg);
    throwAsCaller(MException(eid, msg))
end

end
