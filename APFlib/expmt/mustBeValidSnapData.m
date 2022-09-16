function mustBeValidSnapData(sdata)
%MUSTBEVALIDSNAPDATA Validates a set of snapshot data
%   MUSTBEVALIDSNAPDATA(sdata)
%   checks the validity of a given snapshot data struct using `checkSnapData`.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of Anticipatory Power Flow, and is covered by the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~checkSnapData(sdata, "beVerbose", false)
    eid = "mustBeValidSnapData:notValidSnapData";
    msg = "Snapshot data must be valid.\n";
    msg = sprintf("%sPlease use `checkSnapData()` to see the errors.\n", msg);
    throwAsCaller(MException(eid, msg))
end

end
