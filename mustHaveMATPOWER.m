function mustHaveMATPOWER(~)
%MUSTHAVEMATPOWER Validates MATPOWER installation
%   MUSTHAVEMATPOWER()
%   checks if MATPOWER is installed via looking for `mpver` in the MATLAB search path.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of Anticipatory Power Flow, and is covered by the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~exist('mpver', 'file')
    eid = "mustHaveMATPOWER:noMATPOWER";
    msg = "MATPOWER must be installed first.\n";
    throwAsCaller(MException(eid, msg));
end

end
