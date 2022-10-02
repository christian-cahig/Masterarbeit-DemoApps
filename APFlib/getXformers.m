function [xfm, num] = getXformers(mpc)
%GETXFORMERS Returns the indices of transformers
%   [xfm, num] = GETXFORMERS(mpc)
%   returns the indices of branches representing power transformers, and optionally, the number of
%   power transformers, given a MATPOWER case struct assumed to be using internal indexing.
%
%   Notes
%   ---------------------
%   1. This function relies heavily on the MATPOWER package. It has been developed and tested on
%      version 7.1. The accompanying User's Manual is freely available at
%      https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%      format and indexing, please refer to Appendix B and Section 9.4.
%   2. A '(power) transformer' is defined here as a branch whose corresponding column-9 entry of
%      the branch data matrix, `mpc.branch`, is nonzero.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc {mustBeInterIdxMPC}
end

%% Main
xfm = find(mpc.branch(:, 9) ~= 0);
if nargout > 1, num = length(xfm); end

end
