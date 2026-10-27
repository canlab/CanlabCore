function M = apply_diagonal_correction(M, group_idx, fold_idx, mode)
% apply_diagonal_correction  Replace diagonal entries to remove within-fold bias.
%
% Reproduces (in generalized form) the patch from `generate_RSA3.m`: replace
% each diagonal cell with the mean of off-diagonal same-group cells, so that
% downstream cell-extraction (rsm.cells, rsm.contrast) does not pull a self-
% similarity (= 1 for correlation RSMs) into the same-group bin.
%
% Note on semantics
% -----------------
% The original `generate_RSA3.m` recipe additionally constrained the
% averaged cells to come from different sessions (`ses_numbers(r1) ~=
% ses_numbers(r2)`). That filter relies on per-image session metadata
% which is lost after condition aggregation. In practice, since
% off-diagonal cells already compare DIFFERENT conditions (and therefore
% different images), the "same group, off-diagonal" mean is the natural
% generalization and matches the workflow's downstream cell-extraction
% semantics. The cross-session constraint is only meaningful at the
% image-level RSM, which `compute_rsm` exposes via level='image'.
%
% Inputs
% ------
%   M          [k x k] or [k x k x N] similarity matrix
%
%   group_idx  [k x 1] integer condition labels (1..n_groups). Each entry i
%              identifies which "diagonal-correction group" row/col i belongs
%              to. Diagonals will be corrected per group.
%
%   fold_idx   [k x 1] integer fold labels (1..n_folds) -- IGNORED in the
%              'same_group_offdiag_mean' / 'across_session_mean' modes since
%              fold information is not available at the condition level.
%              Retained for API compatibility and for the image-level path
%              when level='image' is used. Pass NaN(k,1) if unused.
%
%   mode       'same_group_offdiag_mean' (default) -- mean of off-diagonal
%              same-group cells. RECOMMENDED. Use this for condition-aggregated
%              RSMs.
%
%              'across_session_mean' -- alias for 'same_group_offdiag_mean'.
%              Kept for backwards compatibility with the original workflow
%              naming.
%
%              'image_level_across_session' -- requires fold_idx; replaces
%              each diagonal with the mean of cells (i,j) where group(i) ==
%              group(j) and fold(i) != fold(j). Only meaningful at image-level
%              RSM (i.e. when compute_rsm was called with level='image').
%
%              'nan' -- set each same-group diagonal entry to NaN. Useful when
%              the user wants downstream code to apply its own handling.
%
% Output
% ------
%   M          corrected matrix, same size as input

if nargin < 4 || isempty(mode), mode = 'same_group_offdiag_mean'; end
mode = lower(char(mode));
if strcmp(mode, 'across_session_mean'), mode = 'same_group_offdiag_mean'; end

[k1, k2, N] = size(M);
assert(k1 == k2, 'apply_diagonal_correction:M must be square in first two dims.');
k = k1;
group_idx = group_idx(:);
assert(numel(group_idx) == k, ...
    'apply_diagonal_correction: group_idx must be length k = %d.', k);

if nargin < 3 || isempty(fold_idx), fold_idx = nan(k, 1); end
fold_idx = fold_idx(:);

groups = unique(group_idx(~isnan(group_idx)));

for n = 1:N
    slice = M(:, :, n);

    for g_i = 1:numel(groups)
        idx = find(group_idx == groups(g_i));
        if numel(idx) < 2, continue; end   % singletons can't be corrected

        switch mode
            case 'same_group_offdiag_mean'
                % Pull off-diagonal cells within this group
                sub = slice(idx, idx);
                off_diag_mask = ~eye(numel(idx), 'logical');
                vals = sub(off_diag_mask);
                vals = vals(isfinite(vals));
                if isempty(vals), new_val = NaN; else, new_val = mean(vals); end

            case 'image_level_across_session'
                % Pull only cells where folds differ
                vals = [];
                for a = 1:numel(idx)
                    for b = 1:numel(idx)
                        if a == b, continue; end
                        if isnan(fold_idx(idx(a))) || isnan(fold_idx(idx(b)))
                            continue
                        end
                        if fold_idx(idx(a)) ~= fold_idx(idx(b))
                            vals(end+1) = slice(idx(a), idx(b)); %#ok<AGROW>
                        end
                    end
                end
                vals = vals(isfinite(vals));
                if isempty(vals), new_val = NaN; else, new_val = mean(vals); end

            case 'nan'
                new_val = NaN;

            otherwise
                error('apply_diagonal_correction:badMode', ...
                    ['mode must be one of {same_group_offdiag_mean (default), ', ...
                     'across_session_mean (alias), image_level_across_session, nan}; got %s.'], mode);
        end

        for r = 1:numel(idx)
            slice(idx(r), idx(r)) = new_val;
        end
    end

    M(:, :, n) = slice;
end

end
