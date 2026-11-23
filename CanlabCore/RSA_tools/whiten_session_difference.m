function Xw_folds = whiten_session_difference(Xfolds)
% whiten_session_difference  Diagonal whitening using session-difference noise.
%
% Generalizes the recipe from `generate_RSA_accept_crossnobis.m`:
%   noise = X1 - X2          (per voxel, the residual across the two folds)
%   voxel_var = var(noise)   (1 x voxels)
%   X_f_w = X_f ./ sqrt(voxel_var)
%
% Generalized to k folds: noise estimated as the residual of each fold from
% the across-fold mean, then voxel-wise variance pooled across folds.
%
% Input
%   Xfolds    {n_folds x 1} cell, each cell is [k x voxels]
%
% Output
%   Xw_folds  {n_folds x 1} cell, each cell is [k x voxels], whitened

n_folds = numel(Xfolds);
[~, n_vox] = size(Xfolds{1});

% Across-fold mean per condition
M = mean(cat(3, Xfolds{:}), 3, 'omitnan');   % [k x voxels]

% Pool residuals across folds
resid_stack = zeros(0, n_vox);
for f = 1:n_folds
    R = Xfolds{f} - M;
    resid_stack = [resid_stack; R]; %#ok<AGROW>
end

voxel_var = var(resid_stack, 0, 1, 'omitnan');
voxel_var(~isfinite(voxel_var) | voxel_var <= eps) = 1;

Xw_folds = cell(n_folds, 1);
for f = 1:n_folds
    Xw_folds{f} = Xfolds{f} ./ sqrt(voxel_var);
end

end
