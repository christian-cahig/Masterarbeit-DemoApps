function mustBeEqualSize(a, b)
%MUSTBEEQUALSIZE Validates if two arrays have equal sizes
%   MUSTBEEQUALSIZE(a, b)
%   checks if two arrays `a` and `b` have identical sizes. The built-in functions `isequal` and
%   `size` are used to perform the check.
%
%   This is a minimally-changed version of the eponymous custom validation function used in
%   subsection 'Specific Restrictions Using Validation Functions' of the documentation article
%   "Function Argument Validation". The article, part of the MATLAB documentation, is available at
%   https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html.
%
%   Anticipatory Power Flow
%   Copyright (C) 2022 - present, Christian Cahig
%   https://doi.org/10.5281/zenodo.7077324
%
%   This file is part of Anticipatory Power Flow, and is covered by the CC-BY-4.0 License.
%   See `LICENSE` file in the source directory for details.

if ~isequal(size(a),size(b))
    eid = "mustBeEqualSize:notEqualSize";
    msg = "Size of first array must equal that of the second.";
    throwAsCaller(MException(eid, msg))
end

end
