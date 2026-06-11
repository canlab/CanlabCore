function [vw_nuisance, vw_nuisance_comps] = canlab_extract_ventricle_wm_timeseries(mask_image_dir, imgs, varargin)
% canlab_extract_ventricle_wm_timeseries Extract ventricle and white-matter nuisance time series for fMRI.
%
% :Usage:
% ::
%
%     [vw_nuisance, vw_nuisance_comps] = canlab_extract_ventricle_wm_timeseries(mask_image_dir, imgs, ['noplot'])
%
% Extracts data from ventricles and white matter for each time point in
% a functional time series. You can include vw_nuisance in your
% nuisance matrix in first-level model estimation, and/or use the PCA
% components in vw_nuisance_comps (e.g., as in aCompCor).
%
% 5/4/2012 by Tor Wager and Wani Woo.
%
% :Inputs:
%
%   **mask_image_dir:**
%        Directory containing mask files for ventricle and white
%        matter. You can use canlab_create_wm_ventricle_masks.m to
%        produce these mask files. The mask file names must be
%        ventricles.img and white_matter.img. Example:
%        mask_image_dir = fullfile(subjdir, 'Structural/SPGR');
%
%   **imgs:**
%        Image file names that you want to extract data from. Example:
%        imgs = filenames(fullfile(subjdir, ...
%        'Functional/Preprocessed/run1/wra*nii'), 'absolute', 'char');
%
% :Optional Inputs:
%
%   **'noplot':**
%        Suppress the diagnostic plot. The default behavior is to plot
%        the output variables.
%
% :Outputs:
%
%   **vw_nuisance:**
%        T x 2 matrix of MEAN ventricle and white-matter signal.
%        Column 1: ventricle signal. Column 2: white-matter signal.
%        Rows correspond to time points in the input time series.
%
%   **vw_nuisance_comps:**
%        T x 10 matrix containing the first 5 PCA components of the
%        ventricle voxel time series followed by the first 5 PCA
%        components of the white-matter voxel time series.
%
% :Examples:
% ::
%
%     mask_image_dir = fullfile(subjdir, 'Structural/SPGR');
%     imgs = filenames(fullfile(subjdir, ...
%         'Functional/Preprocessed/run1/wra*nii'), 'absolute', 'char');
%
%     % With diagnostic plot
%     [vw_nuisance, vw_nuisance_comps] = ...
%         canlab_extract_ventricle_wm_timeseries(mask_image_dir, imgs);
%
%     % Without plot
%     [vw_nuisance, vw_nuisance_comps] = ...
%         canlab_extract_ventricle_wm_timeseries(mask_image_dir, imgs, 'noplot');
%
% :See also:
%   - canlab_create_wm_ventricle_masks
%   - fmri_data
%   - filenames
%   - pca

%% --------------------------------------------
% Extract data from ventricles
% ---------------------------------------------

doplot = true;

wh = strcmp(varargin, 'noplot');
if any(wh), doplot = false; varargin(wh) = []; end

vmask = filenames(fullfile(mask_image_dir, 'ventricles.img'), 'char', 'absolute');
wmask = filenames(fullfile(mask_image_dir, 'white_matter.img'), 'char', 'absolute');

vdat = fmri_data(imgs, vmask);

vw_nuisance(:, 1) = mean(vdat.dat, 1)'; % avg across voxels

% save the first 5 components

[~, nuisance_comps] = pca(vdat.dat', 'econ');

vw_nuisance_comps = nuisance_comps(:, 1:5);


% ---------------------------------------------
% Extract data from white matter
% ---------------------------------------------

wdat = fmri_data(imgs, wmask);

vw_nuisance(:, 2) = mean(wdat.dat, 1)'; % avg across voxels

% save the first 5 components

[~, n2] = pca(wdat.dat', 'econ');

vw_nuisance_comps = [vw_nuisance_comps n2(:, 1:5)];


%% plot

if doplot
    
    create_figure('vw_nuisance', 3, 1)
    
    plot(zscore(vw_nuisance));
    title('averages (z-scores)')
    
    subplot(3, 1, 2)
    
    plot(vw_nuisance_comps(:, 1:5))
    title('ventricle components')
    
    subplot(3, 1, 3)
    
    plot(vw_nuisance_comps(:, 6:10))
    title('white matter components')
end

return