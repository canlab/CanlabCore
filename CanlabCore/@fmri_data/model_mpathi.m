function stats = model_mpathi(obj,source_mask,target_mask,varargin)
% Models functional coupling between two brain regions using cross-validated 
% partial least squares (PLS) regression (Original code:model_brain_pathway.m)
%
% ------------------------------------------------------------------------------
% OVERVIEW
% ------------------------------------------------------------------------------
% This script is a streamlined and up-to-date version of the original
% [model_brain_pathway.m] function. Unlike the original implementation, which
% compared traditional functional connectivity (based on ROI-averaged
% signals) with the multivariate Pathway Interaction (MPathI) framework
% described in Kragel et al. (2021, Neuron), the present script focuses
% exclusively on the MPathI approach. In addition, this implementation
% estimates a single directed pathway between one source region and one
% target region, rather than modeling four on-target and off-target
% pathways as in model_brain_pathway.
% 
% This function estimates multivariate pathway-level connectivity between a
% source region (X) and a target region (Y) using Partial Least Squares (PLS).
% The model identifies latent population activity in each region whose
% time series covary maximally across observations (e.g., trials or timepoints).
% 
%
% ------------------------------------------------------------------------------
% obj :
%   fmri_data object containing images (e.g., trials or timepoints)
%
% source_mask :
%   fmri_data object with binary mask defining the source region (X)
%
% target_mask :
%   fmri_data object with binary mask defining the target region (Y)
%
% ------------------------------------------------------------------------------
% OPTIONAL NAME–VALUE PAIRS
% ------------------------------------------------------------------------------
% 'Indices' :
%   Integer vector (n_images × 1) defining cross-validation folds and
%   bootstrap blocks (e.g., subject ID). Default: 10-fold CV.
%
% 'Align' :
%   Perform hyperalignment across subjects (requires 'Indices').
%
% 'nboot' :
%   Number of block bootstrap samples for voxelwise inference on weights.
%
% 'plot' :
%   Plot cross-validated latent correlations.
%
% 'noroi' :
%   Do not return masked ROI data in the output (reduces output size).
%
% ------------------------------------------------------------------------------
% Convention used throughout:
%   X, Y        : [images × voxels]
%   T, U        : [images × 1] latent time series
%   Z, V        : [voxels × 1] spatial weights
% ------------------------------------------------------------------------------
% OUTPUT
% ------------------------------------------------------------------------------
% stats :
%   Structure containing:
%
%   Cross-validated results:
%     • stats.xval.latent_correlations   – corr(T̂, Û) per fold
%     • stats.xval.T_latent_timeseries  – predicted source latent time series
%     • stats.xval.U_latent_timeseries  – predicted target latent time series
%
%   Summary statistics:
%     • stats.overall_xval_r             – correlation across all held-out data
%     • stats.overall_xval_dot           – dot product of latent time series
%
%   Whole-sample estimates:
%     • stats.whole.T_latent_timeseries  - source latent time series
%     • stats.whole.U_latent_timeseries  - target latent time series
%     • stats.whole.Z_weights            – source voxel weights
%     • stats.whole.V_weights            – target voxel weights
%
%   Bootstrap results (optional):
%     • stats.PLS_bootstrap_stats_Z
%     • stats.PLS_bootstrap_stats_V
%
% ------------------------------------------------------------------------------
% INTERPRETING RESULTS
% ------------------------------------------------------------------------------
% • The primary measure of pathway strength is stats.overall_xval_r.
% • Cross-validated correlations reflect out-of-sample predictive coupling.
% • Whole-sample weights should be interpreted descriptively unless
%   supported by bootstrap inference.
%
% ------------------------------------------------------------------------------
% Author and copyright information:
% ------------------------------------------------------------------------------
%
% Copyright (C) 2026 Byeol Kim Lux
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Get defaults and initialize user inputs
if any(strcmp(varargin, 'plot'))
    do_plot = true;
else
    do_plot = false;
end

if any(strcmp(varargin,'Align'))
    do_alignment = true;
else
    do_alignment = false;
end

if any(strcmp(varargin,'nboot'))
    do_boot = true;
    nboot  = varargin{find(strcmp(varargin,'nboot'))+1}; % Get number of bootstrap iterations
else
    do_boot = false; % Default: no bootstrapping
end

if any(strcmp(varargin,'noroi'))
    do_roi = false; 
else
    do_roi = true; % Default: save masked ROI data
end

% Default 10-fold cross-validation unless custom indices are provided
% Custom holdout set is specified by a vector of integers, 1 per holdout set
indices=crossvalind('Kfold',size(obj.dat,2),10); 
if any(strcmp(varargin,'Indices'))
    indices = varargin{find(strcmp(varargin,'Indices'))+1}; % Use user-provided indices, like run-folds
end

flip_maps=true; % Flip map orientation (default)

ndim=1; %number of latent dims for PLS

%% convert data and masks to masked objects
% Apply masks to create masked data objects for source and target regions
source_obj=apply_mask(obj,source_mask);
target_obj=apply_mask(obj,target_mask);

%% Cross-validated estimation of latent pathway connectivity
% For each fold:
%   1) Fit PLS on training data
%   2) Derive voxelwise weights (Z, V)
%   3) Predict latent time series in held-out data
%   4) Compute corr(T_hat, U_hat)

% === Prepare Data for Pattern-Based Connectivity Analysis ===

% Convert data into cell arrays for cross-validation
% Each cell contains a Voxels x Images matrix for the corresponding fold
% For 10-fold, images will be devided into 10 cells


for xval_i=1:max(indices)
    un_inds=unique(indices);

    to_align_dat_source{xval_i}=source_obj.dat(:,indices==xval_i);
    to_align_dat_target{xval_i}=target_obj.dat(:,indices==xval_i);
    
    if do_alignment
        % Call Manning hyperalign.m
        % (Currently hyperaligns rows)
        aligned_dat_source = hyperalign(to_align_dat_source{un_inds~=xval_i});
        aligned_dat_target = hyperalign(to_align_dat_target{un_inds~=xval_i});
        
        % X and Y for training [images x voxels] 
        X_train=[aligned_dat_source{:}]';
        Y_train=[aligned_dat_target{:}]';
        % X and Y for test  [images x voxels]
        [~, X_test] = procrustes(mean(cat(3,to_align_dat_source{un_inds~=xval_i}),3), to_align_dat_source{un_inds==xval_i});
        [~, Y_test] = procrustes(mean(cat(3,to_align_dat_target{un_inds~=xval_i}),3), to_align_dat_target{un_inds==xval_i});
        X_test = X_test';
        Y_test = Y_test';
    else
        % X and Y for training [images x voxels]
        X_train=[to_align_dat_source{un_inds~=xval_i}]';
        Y_train=[to_align_dat_target{un_inds~=xval_i}]';
        % X and Y for test [images x voxels]
        X_test=[to_align_dat_source{un_inds==xval_i}]';
        Y_test=[to_align_dat_target{un_inds==xval_i}]';
    end
    

    % === Running PLS regression on the training dataset ===
    % The PLS models estimate the correlation between latent variables (T and U).
    % Xscores (T) and Yscores (U) from plsregress correspond to the latent time series of source and target regions.
    % T, xs, latent time series of the source region (X) [time x 1].
    % U, ys, latent time series of the target region (Y) [time x 1].
    [~,~,T_xval, U_xval] = plsregress(X_train,Y_train,ndim);

    % === Calculate spatial patterns ===
    % Spatial patterns defined as U = X_train * Z', T = Y_train * V' in the
    % Neuron paper. In this code, Z and V are transposed, but no problem.
    % Therefore, Z = pinv(X) * U, V = pinv(Y) * T. 
    % Z represents the spatial pattern for the source (X) [voxel# of X x 1].
    % V represents the spatial pattern for the target (Y) [voxel# of Y x 1].
    Z_xval = pinv(X_train) * U_xval; 
    V_xval = pinv(Y_train) * T_xval;
    
    % === Normalize the test data ===
    % Center each column of the test datasets by subtracting its mean.
    for i=1:size(Y_test,2)
        Y_test(:,i)=Y_test(:,i)-mean(Y_test(:,i));
    end
    
    for i=1:size(X_test,2)
        X_test(:,i)=X_test(:,i)-mean(X_test(:,i));
    end
        
     
    % === Calculate predicted latent time series (T_hat and U_hat) ===
    % T_hat: Predicted latent time series of the source region (X) covarying with the target region (Y).
    % Linear combination of target test data (Y) and spatial patterns of X (V)
    T_hat_xval = Y_test * V_xval;

    % U_hat: Predicted latent time series of the target region (Y) covarying with the source region (X).
    % Linear combination of source test data (X) and spatial patterns of Y (Z)
    U_hat_xval = X_test * Z_xval; 

    % Store predicted latent time series for each region
    stats.xval.T_latent_timeseries(indices==xval_i, 1) = T_hat_xval;
    stats.xval.U_latent_timeseries(indices==xval_i, 1) = U_hat_xval;
    

    % === Calculate latent correlation: Corr(T_hat, U_hat) ===

    % Correlation btw T hat and U hat based on the optimized pair
    latent_correlation(xval_i,:)=corr(T_hat_xval, U_hat_xval);
      
end % End of Cross-validation loop 

% === Cross-validation results ===
stats.xval.whfolds = indices;
stats.xval.latent_correlations = latent_correlation';

%% Model parameters based on whole sample

% Perform hyperalignment across subjects

if do_alignment
    
    % Convert data into cell arrays, with one cell per subject
    % Each cell contains a Voxels x Time (Images) matrix
    for xval_i=1:max(indices)
        to_align_dat_source{xval_i}=source_obj.dat(:,indices==xval_i);
        to_align_dat_target{xval_i}=target_obj.dat(:,indices==xval_i);
    end
    
    % Call Manning hyperalign.m
    % (Currently hyperaligns rows)
    aligned_dat_source = hyperalign(to_align_dat_source{:});
    aligned_dat_target = hyperalign(to_align_dat_target{:});
    
    for xval_i=1:max(indices)
        source_obj.dat(:,indices==xval_i)=aligned_dat_source{xval_i};
        target_obj.dat(:,indices==xval_i)=aligned_dat_target{xval_i};
    end
end

% Fit an overall model for each pathway using the complete dataset, including both training and test data.
% T, latent time series of the source region (X) [time x 1].
% U, latent time series of the target region (Y) [time x 1].
% xl: loading/coefficient for the source region (X), calculated for sign checking.
[XL, ~, T_whole, U_whole] = plsregress(source_obj.dat',target_obj.dat', ndim);

% Flip the latent variable signs to ensure consistency
% Determine the sign flip based on spatial correlation between latent X loadings and source region averages
if flip_maps
    % Enforce consistent sign across folds and whole-sample model.
    % Latent variables are sign-indeterminate in PLS; we align them to
    % positively correlate with the mean source-region signal.
    sign_pathway = sign(corr(mean(XL, 2), mean(source_obj.dat, 2)));
    
    % Adjust scores for X (source) and Y (target) latent variables based on sign flips
    U_whole = U_whole * sign_pathway;
    T_whole = T_whole * sign_pathway;
  
    % Adjust overall latent time series to match the sign flip
    stats.xval.T_latent_timeseries = stats.xval.T_latent_timeseries * sign_pathway;
    stats.xval.U_latent_timeseries = stats.xval.U_latent_timeseries * sign_pathway;
    stats.pathway_sign = sign_pathway;
end

stats.overall_xval_r = corr(stats.xval.T_latent_timeseries, stats.xval.U_latent_timeseries);
stats.overall_xval_dot = dot(stats.xval.T_latent_timeseries, stats.xval.U_latent_timeseries);

stats.whole.T_latent_timeseries = T_whole;
stats.whole.U_latent_timeseries = U_whole;

% Compute latent weights (Z) for source regions and pattern weights (V) for
% target regions with updated T and U.
% Reminder: Z = pinv(X) * U, V = pinv(Y) * T [Voxels x 1]
stats.whole.Z_weights = pinv(source_obj.dat') * U_whole; 
stats.whole.V_weights = pinv(target_obj.dat') * T_whole; 


% Bootstrap inference:
% Blocks defined by 'Indices' are resampled with replacement.
% This preserves within-subject dependence while estimating
% variability of voxelwise pathway weights.
if do_boot    % Perform bootstrap resampling if specified
    bs_source_dat=source_obj.dat';
    bs_target_dat=target_obj.dat';
    for i=1:nboot
        % Resample indices with replacement at the fold level
        [~,rand_subs]=datasample(1:max(indices), max(indices)); %randomly replace whole blocks with replacement
        count_vec=1:length(indices);

        rand_inds=[];
        for ii=1:max(indices)
            rand_inds=[rand_inds count_vec(indices==rand_subs(ii))];
        end

        % Compute bootstrap estimates for source and target weights
        bs_V(i,:)= bootPLS_target_pattern_weights(bs_source_dat(rand_inds,:),bs_target_dat(rand_inds,:),ndim,flip_maps);
        bs_Z(i,:)= bootPLS_source_pattern_weights(bs_source_dat(rand_inds,:),bs_target_dat(rand_inds,:),ndim,flip_maps);
    end
    
   
    % Compute bootstrap statistics for Z and V (weights)
    bs_Z=(mean(bs_Z)./std(bs_Z))';
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    bs_stat=statistic_image;
    bs_stat.volInfo=source_obj.volInfo;
    bs_stat.dat=stats.whole.Z_weights;
    bs_stat.p=bs_P;
    bs_stat.removed_voxels=source_obj.removed_voxels;
    stats.PLS_bootstrap_stats_Z(1)=bs_stat;
        
    bs_Z=(mean(bs_V)./std(bs_V))';
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    bs_stat=statistic_image;
    bs_stat.volInfo=target_obj.volInfo;
    bs_stat.dat=stats.whole.V_weights;
    bs_stat.p=bs_P;
    bs_stat.removed_voxels=target_obj.removed_voxels;
    stats.PLS_bootstrap_stats_V(1)=bs_stat;
end


%% Output data objects
if do_roi
    stats.source_obj=source_obj;
    stats.target_obj=target_obj;
end

%% Optional plot
if do_plot
    % Plot the correlations among PLS-optimized latent timeseries:
    create_figure('MPathI correlations');

    barplot_columns(stats.xval.latent_correlations', 'colors',[1 .7 .3] , 'nofigure');
    title('Cross-validation correlations')
    set(gcf,'color','w','position', [500 500 262 313]);
    set(gca,'fontsize',15);
end


end % main function



%% subfunctions for finding patterns across voxels associated with latent scores in each brain region
function V = bootPLS_target_pattern_weights(X,Y,ndim,flip_maps) %voxels in target region that predict latent variables in source
[xl,~,xs] = plsregress(X,Y,ndim);

if flip_maps
    xs=xs*sign(corr(mean(xl,2),mean(X)'));
end

V = pinv(Y) * xs;
end


function Z = bootPLS_source_pattern_weights(X,Y,ndim,flip_maps) %voxels in source region that predict latent variables in target

[xl,~,~,ys] = plsregress(X,Y,ndim);
if flip_maps
    ys=ys*sign(corr(mean(xl,2),mean(X)'));
end

Z = pinv(X) * ys;

end

