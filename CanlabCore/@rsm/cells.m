function vals = cells(obj, a, b, varargin)
% cells  Extract a vector of cell values from an rsm.
%
% Given two row/column groupings A and B, pulls the (i, j) cells of the
% rsm where i is in A and j is in B. Returns a per-replicate vector of
% reduced (mean / median / sum) values — one scalar per slice of the
% replicate axis (3rd dim).
%
% This is the Phase 2 "down payment" that eliminates the most common
% manual-loop footgun: `mean(M(rows, cols))` returns the COLUMN means of
% the submatrix (a row vector), NOT a scalar. cells() returns a clean
% per-replicate scalar so downstream `ttest(within, between)` calls work.
%
% Within vs between semantics
% ---------------------------
%   A == B   (same grouping for rows and columns):
%       Within-group. Pulls the upper triangle of the square submatrix,
%       excluding the diagonal. n_cells = numel(A) * (numel(A) - 1) / 2.
%
%   A != B   (different groupings):
%       Between-group. Pulls the full rectangular submatrix.
%       n_cells = numel(A) * numel(B).
%
% Grouping inputs
% ---------------
% Each of a, b can be:
%   - char/string -- looked up in obj.groupings
%   - numeric vector -- 1-based row/col indices
%   - logical vector -- length k mask
%
% Usage examples
% --------------
%   % Within-condition
%   v = R.cells('hot', 'hot');                    % [N_subjects x 1]
%
%   % Between-condition
%   v = R.cells('hot', 'warm');                   % [N_subjects x 1]
%
%   % With explicit indices
%   v = R.cells(1:8, 9:16);
%
%   % Suppress Fisher-z (default is on for correlation-like metrics)
%   v = R.cells('hot', 'hot', 'transform','none');
%
%   % Median reduction (default is mean)
%   v = R.cells('hot', 'hot', 'reduction','median');
%
%   % Return ALL cells per replicate instead of reducing
%   V = R.cells('hot', 'warm', 'reduce', false);  % [N_cells x N_subjects]
%
% Optional args
% -------------
%   'transform'   'auto' (default) -- atanh for correlation/spearman/cosine RSMs,
%                                     identity for everything else.
%                 'fisherz'         -- force atanh
%                 'none'            -- raw values
%   'reduction'   'mean' (default) | 'median' | 'sum'
%   'reduce'      logical (default true). If false, returns the full cell
%                 matrix [N_cells x N_replicates] instead of reducing per
%                 replicate.
%
% Output
% ------
%   vals  [N_replicates x 1] when reduce=true (default)
%         [N_cells x N_replicates] when reduce=false

if builtin('numel', obj) > 1
    error('rsm:cells:nonScalar', 'cells() expects a scalar rsm.');
end

p = inputParser;
p.addParameter('transform', 'auto', @(x) ischar(x) || isstring(x));
p.addParameter('reduction', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('reduce',    true,   @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

idx_a = resolve_grouping(obj, a, 'a');
idx_b = resolve_grouping(obj, b, 'b');

k = size(obj.dat, 1);
N = size(obj.dat, 3);

is_within = isequal(sort(idx_a(:)), sort(idx_b(:)));

% Build the cell mask once
if is_within
    sub_mask = triu(true(numel(idx_a)), 1);     % exclude diagonal
else
    sub_mask = true(numel(idx_a), numel(idx_b));
end
n_cells = nnz(sub_mask);

if n_cells == 0
    warning('rsm:cells:empty', 'cells() returned no cells (singleton within-grouping?).');
    if opt.reduce, vals = nan(N, 1); else, vals = nan(0, N); end
    return
end

% Decide transform
transform = lower(char(opt.transform));
if strcmp(transform, 'auto')
    if obj.is_dissimilarity || ~ismember(lower(obj.metric), {'correlation','spearman','cosine'})
        transform = 'none';
    else
        transform = 'fisherz';
    end
end

% Pull cells per replicate
all_cells = zeros(n_cells, N);
for n = 1:N
    slice = obj.dat(:, :, n);
    sub = slice(idx_a, idx_b);
    all_cells(:, n) = sub(sub_mask);
end

% Apply transform
switch transform
    case 'fisherz'
        all_cells(all_cells >  0.9999999) =  0.9999999;
        all_cells(all_cells < -0.9999999) = -0.9999999;
        all_cells = atanh(all_cells);
    case 'none'
        % no-op
    otherwise
        error('rsm:cells:badTransform', ...
            'transform must be ''auto'', ''fisherz'', or ''none''; got %s.', transform);
end

if ~logical(opt.reduce)
    vals = all_cells;
    return
end

switch lower(char(opt.reduction))
    case 'mean',   vals = mean(all_cells, 1, 'omitnan')';
    case 'median', vals = median(all_cells, 1, 'omitnan')';
    case 'sum',    vals = sum(all_cells, 1, 'omitnan')';
    otherwise
        error('rsm:cells:badReduction', ...
            'reduction must be ''mean'', ''median'', or ''sum''; got %s.', opt.reduction);
end

end


function idx = resolve_grouping(obj, x, label)
% Coerce a grouping spec (name / numeric / logical) to a 1-based index vector.

k = size(obj.dat, 1);

if ischar(x) || isstring(x)
    name = char(x);
    if ~isfield(obj.groupings, name)
        avail = fieldnames(obj.groupings);
        if isempty(avail)
            error('rsm:cells:noGroupings', ...
                ['No .groupings defined on this rsm; cannot resolve grouping name "%s". ', ...
                 'Pass numeric indices or set obj.groupings before calling cells().'], name);
        end
        error('rsm:cells:unknownGrouping', ...
            'Unknown grouping name "%s" for argument %s. Available: %s', ...
            name, label, strjoin(avail, ', '));
    end
    idx = obj.groupings.(name);
elseif islogical(x)
    if numel(x) ~= k
        error('rsm:cells:badLogical', ...
            'Logical grouping %s must be length k=%d; got %d.', label, k, numel(x));
    end
    idx = find(x);
elseif isnumeric(x)
    idx = x;
else
    error('rsm:cells:badGrouping', ...
        'Grouping %s must be a name, numeric index, or logical mask.', label);
end

idx = idx(:);
if any(idx < 1) || any(idx > k)
    error('rsm:cells:outOfRange', ...
        'Grouping %s contains indices outside 1..%d.', label, k);
end

end
