function [vw_nuisance, vw_nuisance_comps] = canlab_extract_ventricle_wm_timeseries(mask_image_dir, imgs, varargin)

% function [nuisance, nuisance_comps] = canlab_extract_ventricle_wm_timeseries(mask_image_dir)
%
% This function extracts data from ventricles and white matter for each
% time series. You can include vw_nuisance in your nuisance matrix in the
% first level model estimation. 
%
% output: "vw_nuisance" returns the MEAN ventricle and white matter signal
%         1st column: ventricle signal, 2nd column: white matter signal, 
%         rows: time series
%
%         "vw_nuisance_comps" returns first 5 components from PCA.
%
% input: 
% mask_image_dir: a directory that has mask files for ventricle and white matter.
%        You can use "canlab_create_wm_ventricle_masks.m" to get these mask files.
%        The mask file names must be "ventricles.img' and 'white_matter.img'. 
%        e.g) mask_image_dir = fullfile(subjdir, 'Structural/SPGR');
%
% imgs: image file names that you want to extract data from. 
%        e.g) imgs = filenames(fullfile(subjdir, 'Functional/Preprocessed/run1/wra*nii'), 'absolute', 'char');
%
% option (varargin): 'noplot' - the default setting is plotting output
%                    variables. You can use 'noplot' option if you don't want the plot. 
%
% 5/4/2012 by Tor Wager and Wani Woo

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

[~, nuisance_comps] = princomp(vdat.dat', 'econ');

vw_nuisance_comps = nuisance_comps(:, 1:5);


% ---------------------------------------------
% Extract data from white matter
% ---------------------------------------------

wdat = fmri_data(imgs, wmask);

vw_nuisance(:, 2) = mean(wdat.dat, 1)'; % avg across voxels

% save the first 5 components

[~, n2] = princomp(wdat.dat', 'econ');

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