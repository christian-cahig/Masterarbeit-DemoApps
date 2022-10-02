function [] = installAPFlib(~)
%INSTALLAPFLIB Installs the APF library
%   INSTALLAPFLIB()
%   adds `APFlib/`, `APFlib/utils/`, and `APFlib/expmt/` to the MATLAB search path.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of a repository under the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

%% Arguments
arguments
    ~
end

%% Utilities
INDENT = "   ";

%% Main
% Check if MATPOWER is available
mustHaveMATPOWER();

% Check if CVX is available
mustHaveCVX();

% Add directories to path
addpath('APFlib', '-end');
addpath('APFlib/utils/', '-end');
addpath('APFlib/expmt/', '-end');

%% Summary
fprintf("%sStuff are ready to use.\n", INDENT);

end