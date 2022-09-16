function [flag] = isMPC(mpc)
%ISMPC Checks if a struct defines a MATPOWER case
%   flag = ISMPC(mpc)
%   returns a logical, indicating if the given struct defines a MATPOWER case.
%
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf.
%   2. For `mpc` to be a MATPOWER case, it must be a struct with the fields `baseMVA`, `bus`,
%      `branch`, `gen`, `gencost`, and `version`. For more details, please refer to Appendix B of
%      the MATPOWER version 7.1 User's Manual.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of Anticipatory Power Flow, and is covered by the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc struct
end

%% Main
flag = all(isfield(mpc, {'baseMVA', 'bus', 'branch', 'gen', 'gencost', 'version'}));

end
