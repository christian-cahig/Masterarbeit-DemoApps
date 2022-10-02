function [sdata] = prepBaseData(mpc)
%PREPBASEDATA Prepares a built-in MATPOWER case data according to Appendix B.2 of the manuscript
% 
%   This function relies heavily on the MATPOWER package. It has been developed and tested on
%   version 7.1. The accompanying User's Manual is freely available at
%   https://matpower.org/docs/MATPOWER-manual-7.1.pdf. For more details on MATPOWER's data file
%   format and indexing, please refer to Appendix B and Section 9.4.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    mpc struct = case57
end

%% Main
sdata = mpc;
sdata.bus(:, 3) = max(0, sdata.bus(:, 3));
sdata.gen(:, 8) = 1;
sdata.branch(:, 11) = 1;
sdata.gencost(:, 1) = 2;
sdata.gencost = sdata.gencost(:, 1:7);
sdata.gencost(:, 5) = 1;
sdata.gencost(:, 6) = 1;
sdata.gencost(:, 7) = 0;

end