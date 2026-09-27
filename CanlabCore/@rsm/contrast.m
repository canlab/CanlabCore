function vals = contrast(obj, varargin)
% contrast  Compute a per-replicate contrast scalar from cell groupings.
%
% A contrast compares one set of cells (A) against another (B): for each
% replicate, it computes (reduction of cells in A) minus (reduction of
% cells in B). For one-sample tests against zero, pass only A.
%
% Usage
% -----
%   % One-sample (cells_A only): mean of A's cells per replicate
%   v = R.contrast('hot');
%   v = R.contrast({'hot','hot'});
%
%   % Two-sample (paired): mean of A - mean of B, per replicate
%   v = R.contrast('hot', 'warm');                              % within-hot vs within-warm
%   v = R.contrast({'hot','warm'}, {'hot','imagine'});          % HW vs HI
%
%   % With explicit indices
%   v = R.contrast(1:8, 9:16);
%
% Each "cells spec" can be:
%   - char/string                 -> within-group (R.cells(name, name))
%   - {a, b} 1x2 cell             -> between-group (R.cells(a, b))
%   - numeric/logical index       -> within-group on that index set
%
% Optional name-value
%   'transform'   'auto' (default) | 'fisherz' | 'none'
%   'reduction'   'mean' (default) | 'median' | 'sum'
%
% Output
%   vals  [N_replicates x 1] column vector

if builtin('numel', obj) > 1
    error('rsm:contrast:nonScalar', 'contrast() expects a scalar rsm.');
end

if isempty(varargin)
    error('rsm:contrast:noArgs', 'Need at least one cells spec.');
end

% Detect: 1 cells-spec (one-sample) vs 2 cells-specs (two-sample)
% Cells specs may be followed by name-value pairs ('transform', ...).
% We assume the first 1 or 2 positional args are specs and the rest are NV pairs.

n_pos = 0;
for i = 1:numel(varargin)
    if ischar(varargin{i}) || isstring(varargin{i})
        s = char(varargin{i});
        if any(strcmpi(s, {'transform','reduction','reduce'}))
            break
        end
    end
    n_pos = n_pos + 1;
    if n_pos == 2, break; end
end

if n_pos == 0
    error('rsm:contrast:noSpec', 'Need at least one cells spec before name-value pairs.');
end

specA = varargin{1};
nv = varargin(n_pos+1:end);

if n_pos == 1
    vA = pull_cells(obj, specA, nv{:});
    vals = vA;
else
    specB = varargin{2};
    vA = pull_cells(obj, specA, nv{:});
    vB = pull_cells(obj, specB, nv{:});
    vals = vA - vB;
end

end


function v = pull_cells(obj, spec, varargin)
% Resolve a cells spec and pull values via R.cells.
if ischar(spec) || isstring(spec)
    a = char(spec); b = a;
elseif iscell(spec) && numel(spec) == 2
    a = spec{1}; b = spec{2};
elseif (isnumeric(spec) || islogical(spec))
    a = spec; b = spec;
else
    error('rsm:contrast:badSpec', ...
        'cells spec must be name, {a,b} cell, or index vector.');
end
v = obj.cells(a, b, varargin{:});
end
