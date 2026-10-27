function Xfolds = build_fold_pattern_matrices(dat_mat, group_idx, fold_idx)
% build_fold_pattern_matrices  Build per-fold [k x voxels] pattern matrices.
%
% Generalizes `build_accept_runfold_matrices` from
% `generate_RSA_accept_crossnobis.m`. For each unique fold f, returns a
% [k x voxels] matrix where row i is the mean pattern for condition i within
% fold f.
%
% Inputs
%   dat_mat    [voxels x n_images] data matrix
%   group_idx  [n_images x 1] integer condition labels (1..k)
%   fold_idx   [n_images x 1] integer fold labels (1..n_folds)
%
% Output
%   Xfolds     {n_folds x 1} cell, each cell is [k x voxels]

n_images = size(dat_mat, 2);
if numel(group_idx) ~= n_images || numel(fold_idx) ~= n_images
    error('build_fold_pattern_matrices:badShape', ...
        'group_idx and fold_idx must each be length size(dat_mat,2) = %d.', n_images);
end

group_idx = group_idx(:);
fold_idx  = fold_idx(:);

folds = unique(fold_idx);
groups = unique(group_idx);
k = numel(groups);
n_folds = numel(folds);
n_vox = size(dat_mat, 1);

Xfolds = cell(n_folds, 1);
for f = 1:n_folds
    X = nan(k, n_vox);
    is_fold = fold_idx == folds(f);
    for g = 1:k
        is_group = group_idx == groups(g);
        rows = is_fold & is_group;
        if any(rows)
            X(g, :) = mean(dat_mat(:, rows), 2, 'omitnan')';
        end
    end
    Xfolds{f} = X;
end

end
