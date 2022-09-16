function mustBeInterIdxMPC(mpc)
%MUSTBEINTERIDXMPC Validates if a MATPOWER case is internally indexed
%   MUSTBEINTERIDXMPC(mpc)
%   checks if a given MATPOWER case struct is adopting MATPOWER's internal indexing. The check is
%   performed via `isInterIdxMPC`.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of Anticipatory Power Flow, and is covered by the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~isInterIdxMPC(mpc)
    eid = "mustBeInterIdxMPC:notInterIdxMPC";
    msg = "Must be an internally-indexed MATPOWER case struct.\n";
    msg = sprintf("%sUse `ext2int()` to convert to internal indexing.\n", msg);
    throwAsCaller(MException(eid, msg))
end

end
