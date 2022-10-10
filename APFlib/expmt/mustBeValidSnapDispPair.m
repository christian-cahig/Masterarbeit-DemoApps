function mustBeValidSnapDispPair(sdata, ddata)
%MUSTBEVALIDSNAPDISPPAIR Validates a pair of snapshot and dispatch data sets
%   MUSTBEVALIDSNAPDISPPAIR(sdata, ddata)
%   checks the validity of a given pair of snapshot and dispatch data structs using the function
%   `checkSnapDispPair`.
%
%   A Computational Approach to Anticipating Supply Injections and Bus Voltages in Steady-State Power System Analysis
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7180586
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~checkSnapDispPair(sdata, ddata, "beVerbose", false)
    eid = "mustBeValidSnapDispPair:notValidSnapDispPair";
    msg = "Pair of snapshot and dispatch data must be valid.\n";
    msg = sprintf("%sPlease use `checkSnapDispPair()` to see the errors.\n", msg);
    throwAsCaller(MException(eid, msg))
end

end
