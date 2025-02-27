function stats = model_brain_pathway(obj,source_one,source_two,target_one,target_two,varargin)
% Models connections between a pair of brain regions using partial least squares regression
%
% Overview:
% This function models the relationship between pairs of brain regions using 
% Partial Least Squares (PLS) regression to identify a set of weights in a source 
% region (X) that explain activity in a target region (Y). 
% This many to many mapping is performed using Partial Least Squares regression, and is based on the
% idea that connected neural populations in the two brain regions produce
% correlated signals (e.g., fMRI activation).
%
%
% Key Features:
% 1. Perform a typical cross-validated connectivity analysis as a comparison condition
%   - Estimate simple correlations using k-fold CV
%   - Analyze four pathways connecting a pair of brain regions, using region averages
%     However, We are only interested in pathways one and four, while pathways two and three will serve as comparison conditions.
%      * Pathway_one:    X1 (Source 1) -> Y1 (Target 1)
%      * Pathway_two:    X2 (Source 2) -> Y1 (Target 1)
%      * Pathway_three:  X1 (Source 1) -> Y2 (Target 2)
%      * Pathway_four:   X2 (Source 2) -> Y2 (Target 2)
%
% 2. Cross-validate PLS model estimating pathways
%   - Use input cross-validation indices if provided
%   - Optionally create hyperalignment model on training data and warp test data 
%     into group mean space
%   - Estimate PLS models for the four pathways with one latent source per pathway.
% 
%   - For each pathway:
%      * Compute Xscores (T) and Yscores (U) from PLS regression on the training data.
%      * Xscores and Yscores represent timeseries of latent population activity covarying between regions.
%      * Calculate Z weights by regressing voxels in X onto latent timeseries of Y (U):
%           Z = pinv(X_train) * U
%      * Calculate V weights by regressing voxels in Y onto latent timeseries of X (T):
%           V = pinv(Y_train) * T
%      * Derive Xscores (T_hat) and Yscores (U_hat) for test data using Z and V weights 
%        from training data:
%           U_hat = X_test * Z, Y_hat = Y_test * V
%      * Compute Pearson's correlation coefficient between latent timeseries. 
%      * Higher correlations are expected for matched pathways (e.g., T_hat and U_hat from X1-Y1) 
%        compared to mismatched pathways (e.g., T_hat from X2-Y1 and U_hat from X1-Y1).
% 
% 3. Estimate the PLS models on the full sample, 
%    - Estimate Xscore, Yscore, V, and Z weights on the full dataset.
%    - Optional bootstrap for significance inference of voxel-level weights (V and Z)
%
% 
% Usage:
%    stats = model_brain_pathway(obj,source_one,source_two,target_one,target_two,varargin);
%
% 
% Default Settings:
% - 10-fold cross-validation (no grouping; see optional inputs) for testing connection strength
% - Block bootstrap with random selection of observations in blocks for significance of pattern weights
% - If running across multiple participants, specify 'Indices' to properly hangle subject grouping
% 
% 
% :Inputs:
%
%   **obj:**
%        An fmri_data image object with one or more images loaded (e.g., single-trial
%        or preprocessed time-series)
%
%   **source_one, source_two:**
%        An fmri_data object whose .dat field is a binary mask specifying
%        voxels that belong to source regions (X1, X2)
%
%   **target_one, target_two:**
%        An fmri_data object whose .dat field is a binary mask specifying
%        voxels that belong to target regions (Y1, Y2)
%
%
% :Optional inputs:
%
%   **'plot'**
%      Create bar plots of output pathway strength
%
%   **'names'**
%      Followed by cell array of pathway names
%      e.g., {'Path1' 'Offtarget1' 'Offtarget2' 'Path2'}
%      e.g., {'L LGN->L V1' 'L LGN->L A1' 'L MGN->L V1' 'L MGN->L A1'};
%
%   **'nboot'**
%      Followed by the number of bootstrap samples to conduct for
%      estimating voxel-wise significance of V and Z weights
%
%   **'Indices'**
%      An integer vector specifying participant (or blocking) ID
%      Used in two ways:
%      (1) to define cross-validation holdout set for each observation, n_obs x 1. 
%       This is useful if you want to, e.g., leave out all images from a participant.
%      (2) to define blocks for block bootstrap
%
%   **'Align'**
%      Perform hyperalignment, requires input 'Indices' to be specified
%      with each index corresponding to a subject
% 
%   **'noroi'**
%      It does not save the masked ROI data from the input, which significantly reduces the size of the output.
% 
% :Outputs:
%
%   **stats:**
%        Structure including:
%           - simple_correlations: correlations results of cross-validated
%           connectivity analysis that used region averages. Pathways are
%           columns 1 and 4.  Off-target pathways are columns 2 and 3.
%           - latent_correlations: correlations among PLS-optimized latent timeseries (T_hat and U_hat)
%           - latent_timeseries: cross-validated time series of the latent scores 
%             latent_timeseries_source: T_hat = Y_test * V
%             latent_timeseries_target: U_hat = X_test * Z
%
%           - Overall cross-validated correlations and dot products for
%           multi-level analysis.  If you run this on a series of
%           individual subjects, each with their own latent pattern, then
%           you can use these as summary statistics for 2nd-level group
%           analyses:
%             stats.path1_overall_xval_r = Correlation between path 1 source and target
%             stats.path2_overall_xval_r = Correlation between path 2 source and target
%             stats.path1_overall_xval_dot = Dot product between path 1 source and target
%             stats.path2_overall_xval_dot = Dot product between path 2 source and target
%
%           - latent_correlation_interaction_ttest: T-test on Fisher-transformed latent_correlation contrasting on-target vs. off-target on two pathways
%             simple_correlation_interaction_ttest: T-test on Fisher-transformed simple_correlation contrasting on-target vs. off-target on two pathways
%             latent_correlation_pathway_*_ttest: T-test on Fisher-transformed latent_correlation contrasting on-target vs. off-target on pathway one or two
%             simple_correlation_pathway_*_ttest: T-test on Fisher-transformed simple_correlation contrasting on-target vs. off-target on pathway one or two
%           - PLS_bootstrap_stats_Z: statistic image for PLS regression coefficients
%           - PLS_bootstrap_stats_V: statistic image for PLS regression coefficients
%           - Z_pathway_one/four: estimated patterns in source regions that covary with their target region
%           - V_pathway_one/four: estimated patterns in target regions that covary with their source region
%           - source_one_obj (optional): resampled and masked fmri data object for source one
%           - source_two_obj (optional): resampled and masked fmri data object for source two
%           - target_one_obj (optional): resampled and masked fmri data object for target one
%           - target_two_obj (optional): resampled and masked fmri data object for target one
%
% :Examples:
% ::
%
% % Example 1: Time series
% % ----------------------------------------------------------
% See canlab_help_MPathI_multivariate_pathway.mlx
% Walkthrough on canlab.github.io
%
%
% % Example 2: Multi-subject toy image series example
% % ----------------------------------------------------------
% % Load multi-subject image object (toy example, 6 images per subject)
% imgs = load_image_set('bmrk3'); 
%
% % Define ROIs
% atlas_obj = load_atlas('canlab2018_2mm');
% vpl = select_atlas_subset(atlas_obj, {'VPL'}, 'flatten');
% pbn = select_atlas_subset(atlas_obj, {'pbn'}, 'flatten');
% cea = select_atlas_subset(atlas_obj, {'Amygdala_CM'}, 'flatten');
% dpins = select_atlas_subset(atlas_obj, {'Ig'}, 'flatten');
%
% % Define participant indices for xval and bootstrap
% indx = imgs.additional_info.subject_id;
%
% % Run it!
% % This estimates PLS models with two pathways: vpl->dpins and pbn->cea
% % It tests the cross-validated connectivity strength in each pathway
% % and across pathways (off=target), and compares target vs off-target pathways
% % 
% stats = model_brain_pathway(imgs, vpl, pbn, dpins, cea, 'Indices', indx);
%
% % Plot the correlations among region averages:
% figure; barplot_columns(stats.simple_correlations, 'names', {'vpl->dpins' 'pbn->dpins' 'vpl->cea' 'pbn->cea'});
%
% % Plot the correlations among PLS-optimized latent timeseries:
% figure; barplot_columns(stats.latent_correlations, 'names', {'vpl->dpins' 'pbn->dpins' 'vpl->cea' 'pbn->cea'});
%
% % Run it again with bootstrapping too:
% stats = model_brain_pathway(imgs, vpl, pbn, dpins, cea, 'Indices', indx, 'nboot', 1000);
%
% :See also: Kragel et al., Neuron 2021.
%

% ..
%    Programmers' notes:
% 10/10/2019 Wrote funcion (Phil Kragel)
% 4/14/2020 Major overhaul and documentation (Phil Kragel and Tor Wager)
% 4/15/2020 Add more documentation (Phil Kragel)
% 9/16/2020 Tor Wager - update documentation
% 10/25/2022 Byeol Lux & Tor Wager - update documentation and add latent_timeseries_pathway
% 12/08/2024 Byeol: Updated some variable names to better align with the Kragel 2021 Neuron paper 
%                   (refer to supplementary figures 1) and enhanced the documentation.
% 2/23/2025 Tor Wager - update documentation and help examples, and 'plot' option
%
% Notes for future development:
% - could use xval_stratified_holdout_leave_whole_subject_out in future versions for cross-validation
% - could write separate hyperalignment object method
% ..
% ..
%     Author and copyright information:
%
%     Copyright (C) 2019 Phil Kragel
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..

%% Get defaults and initialize user inputs

if any(strcmp(varargin, 'names'))
    pathwaynames = varargin{find(strcmp(varargin, 'names'))+1};
    if ~iscell(pathwaynames), error('Enter names of pathways in cell array, {''Path1'' ''Offtarget1'' ''Offtarget2'' ''Path2''}'); end

else
    pathwaynames = {'Path1' 'Offtarget1' 'Offtarget2' 'Path2'};
end
stats.pathwaynames = pathwaynames;

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
source_one_obj=apply_mask(obj,source_one);
source_two_obj=apply_mask(obj,source_two);
target_one_obj=apply_mask(obj,target_one);
target_two_obj=apply_mask(obj,target_two);


%% Perform a typical cross-validated connectivity analysis as a comparison condition

% === Estimate Simple Correlations Using Region Averages ===
% Perform 10-fold cross-validation (default) or user-specified validation
% Compute simple correlations for four pathways connecting two regions:
% Pathways:
% - Pathway 1: Source X1 <-> Target Y1
% - Pathway 2: Source X2 <-> Target Y2


for k=1:max(indices) %for each fold
    
    % Split data into training and testing sets for the current fold
    train = indices~=k; % Training set: all data except current fold
    test = ~train;      % Testing set: data in the current fold
   
    % Simple linear regression between mean region activity for training data
    beta_pathway_one_mean = glmfit(mean(source_one_obj.dat(:,train))',mean(target_one_obj.dat(:,train))');
    % beta_pathway_two_mean = glmfit(mean(source_two_obj.dat(:,train))',mean(target_one_obj.dat(:,train))');  % they have never used
    % beta_pathway_three_mean = glmfit(mean(source_one_obj.dat(:,train))',mean(target_two_obj.dat(:,train))');
    beta_pathway_four_mean = glmfit(mean(source_two_obj.dat(:,train))',mean(target_two_obj.dat(:,train))');

    % Predict target region activity for test data using regression coefficients
    yhat_pathway_one_mean(test,:)=[ones(size(mean(source_one_obj.dat(:,test))',1),1) mean(source_one_obj.dat(:,test))'] * beta_pathway_one_mean; % Y1_hat predicted by X1
    % yhat_pathway_two_mean(test,:)=[ones(size(mean(source_two_obj.dat(:,test))',1),1) mean(source_two_obj.dat(:,test))'] * beta_pathway_two_mean;
    % yhat_pathway_three_mean(test,:)=[ones(size(mean(source_one_obj.dat(:,test))',1),1) mean(source_one_obj.dat(:,test))'] * beta_pathway_three_mean;
    yhat_pathway_four_mean(test,:)=[ones(size(mean(source_two_obj.dat(:,test))',1),1) mean(source_two_obj.dat(:,test))'] * beta_pathway_four_mean; % Y2_hat predicted by X2

    % Calculate simple correlations between predicted and actual target region activity
    % with column order: Y1_hat & Y1, Y1_hat & Y2, Y2_hat & Y1, Y2_hat & Y2
    stats.simple_correlations(k,:)=[corr(yhat_pathway_one_mean(test),mean(target_one_obj.dat(:,test))') corr(yhat_pathway_one_mean(test),mean(target_two_obj.dat(:,test))') ...
        corr(yhat_pathway_four_mean(test),mean(target_one_obj.dat(:,test))') corr(yhat_pathway_four_mean(test),mean(target_two_obj.dat(:,test))')];    

end


%% estimate pattern-based connectivity with kfold or user-specified xval

% === Prepare Data for Pattern-Based Connectivity Analysis ===

% Convert data into cell arrays for cross-validation
% Each cell contains a Voxels x Images matrix for the corresponding fold
% For 10-fold, images will be devided into 10 cells

% note by kodiweera: convert to a cell array before the cross validation loop. I took this loop out of the main loop below to make the the program little faster as this is repeatedlty run during the k-fold validation.
% Convert to cell array, one cell per subject
% Voxels x Images matrix
for kk=1:max(indices)
    to_align_dat_one{kk}=source_one_obj.dat(:,indices==kk);
    to_align_dat_two{kk}=source_two_obj.dat(:,indices==kk);
    to_align_dat_three{kk}=target_one_obj.dat(:,indices==kk);
    to_align_dat_four{kk}=target_two_obj.dat(:,indices==kk);

end


for k=1:max(indices)

    un_inds=unique(indices);
    % for kk=1:max(indices)
    %    to_align_dat_one{kk}=source_one_obj.dat(:,indices==kk);
    %    to_align_dat_two{kk}=source_two_obj.dat(:,indices==kk);
    %    to_align_dat_three{kk}=target_one_obj.dat(:,indices==kk);
    %    to_align_dat_four{kk}=target_two_obj.dat(:,indices==kk);
    % end
    
    if do_alignment
        % Call Manning hyperalign.m
        % (Currently hyperaligns rows)
        aligned_dat_one = hyperalign(to_align_dat_one{un_inds~=k});
        aligned_dat_two = hyperalign(to_align_dat_two{un_inds~=k});
        aligned_dat_three = hyperalign(to_align_dat_three{un_inds~=k});
        aligned_dat_four = hyperalign(to_align_dat_four{un_inds~=k});
        
        % X and Y from train dataset. [voxels x time(images)] 
        source_one_train=[aligned_dat_one{:}];
        source_two_train=[aligned_dat_two{:}];
        target_one_train=[aligned_dat_three{:}];
        target_two_train=[aligned_dat_four{:}];
        
        [~, source_one_test] = procrustes(mean(cat(3,to_align_dat_one{un_inds~=k}),3), to_align_dat_one{un_inds==k});
        [~, source_two_test] = procrustes(mean(cat(3,to_align_dat_two{un_inds~=k}),3), to_align_dat_two{un_inds==k});
        [~, target_one_test] = procrustes(mean(cat(3,to_align_dat_three{un_inds~=k}),3), to_align_dat_three{un_inds==k});
        [~, target_two_test] = procrustes(mean(cat(3,to_align_dat_four{un_inds~=k}),3), to_align_dat_four{un_inds==k});
        
    else
        % X and Y from train dataset. [voxels x time(images)]
        source_one_train=[to_align_dat_one{un_inds~=k}];
        source_two_train=[to_align_dat_two{un_inds~=k}];
        target_one_train=[to_align_dat_three{un_inds~=k}];
        target_two_train=[to_align_dat_four{un_inds~=k}];
        
        source_one_test=[to_align_dat_one{un_inds==k}];
        source_two_test=[to_align_dat_two{un_inds==k}];
        target_one_test=[to_align_dat_three{un_inds==k}];
        target_two_test=[to_align_dat_four{un_inds==k}];
    end
    

    % === Running PLS regression on the training dataset ===
    % The PLS models estimate the correlation between latent variables (T and U).
    % Xscores (T) and Yscores (U) from plsregress correspond to the latent time series of source and target regions.
    % xs: T, latent time series of the source region (X) [time x 1].
    % ys: U, latent time series of the target region (Y) [time x 1].
    [~,~,xs_pathway_one, ys_pathway_one] = plsregress(source_one_train',target_one_train',ndim);
    [~,~,xs_pathway_two, ys_pathway_two] = plsregress(source_two_train',target_one_train',ndim);
    [~,~,xs_pathway_three, ys_pathway_three] = plsregress(source_one_train',target_two_train',ndim);
    [~,~,xs_pathway_four, ys_pathway_four] = plsregress(source_two_train',target_two_train',ndim);

    % === Calculate spatial patterns ===
    % Spatial patterns defined as U = X_train * Z', T = Y_train * V' in the
    % Neuron paper. In this code, Z and V are transposed, but no problem.
    % Therefore, Z = pinv(X) * U, V = pinv(Y) * T. 
    % Z represents the spatial pattern for the source (X) [voxel# of X x 1].
    % V represents the spatial pattern for the target (Y) [voxel# of Y x 1].
    Z_pathway_one = pinv(source_one_train') * ys_pathway_one; 
    V_pathway_one = pinv(target_one_train') * xs_pathway_one;

    Z_pathway_two = pinv(source_two_train') * ys_pathway_two;
    V_pathway_two = pinv(target_one_train') * xs_pathway_two;

    Z_pathway_three = pinv(source_one_train') * ys_pathway_three;
    V_pathway_three = pinv(target_two_train') * xs_pathway_three;

    Z_pathway_four = pinv(source_two_train') * ys_pathway_four;
    V_pathway_four = pinv(target_two_train') * xs_pathway_four;
    
    % === Normalize the test data ===
    % Center each column of the test datasets by subtracting its mean.
    Ytest_target_one=target_one_test';
    for i=1:size(Ytest_target_one,2)
        Ytest_target_one(:,i)=Ytest_target_one(:,i)-mean(Ytest_target_one(:,i));
    end
    
    Ytest_target_two=target_two_test';
    for i=1:size(Ytest_target_two,2)
        Ytest_target_two(:,i)=Ytest_target_two(:,i)-mean(Ytest_target_two(:,i));
    end
    
    Xtest_source_one=source_one_test';
    for i=1:size(Xtest_source_one,2)
        Xtest_source_one(:,i)=Xtest_source_one(:,i)-mean(Xtest_source_one(:,i));
    end
        
    Xtest_source_two=source_two_test';
    for i=1:size(Xtest_source_two,2)
        Xtest_source_two(:,i)=Xtest_source_two(:,i)-mean(Xtest_source_two(:,i));
    end
    
    % optimized models
    % Revision 6/4/2021, Phil Kragel and Tor Wager - documentation
    % Revised by Byeol, changed the variable names based on Kragel 2021 Neuron paper and improve documentation

    % === Calculate predicted latent time series (T_hat and U_hat) ===
    % T_hat: Predicted latent time series of the source region (X) covarying with the target region (Y).
    % Linear combination of target test data (Y) and spatial patterns of X (V)
    T_hat_pathway_one = Ytest_target_one * V_pathway_one;
    T_hat_pathway_two = Ytest_target_one * V_pathway_two;
    T_hat_pathway_three = Ytest_target_two * V_pathway_three;
    T_hat_pathway_four = Ytest_target_two * V_pathway_four;

    % U_hat: Predicted latent time series of the target region (Y) covarying with the source region (X).
    % Linear combination of source test data (X) and spatial patterns of Y (Z)
    U_hat_pathway_one = Xtest_source_one * Z_pathway_one; 
    U_hat_pathway_two = Xtest_source_two * Z_pathway_two;
    U_hat_pathway_three = Xtest_source_one * Z_pathway_three;
    U_hat_pathway_four = Xtest_source_two * Z_pathway_four;

    % Store predicted latent time series for each pathway
    stats.latent_timeseries_source(indices==k, 1) = T_hat_pathway_one(:, 1);
    stats.latent_timeseries_source(indices==k, 2) = T_hat_pathway_two(:, 1);    
    stats.latent_timeseries_source(indices==k, 3) = T_hat_pathway_three(:, 1);
    stats.latent_timeseries_source(indices==k, 4) = T_hat_pathway_four(:, 1);   

    stats.latent_timeseries_target(indices==k, 1) = U_hat_pathway_one(:, 1);
    stats.latent_timeseries_target(indices==k, 2) = U_hat_pathway_two(:, 1);
    stats.latent_timeseries_target(indices==k, 3) = U_hat_pathway_three(:, 1);
    stats.latent_timeseries_target(indices==k, 4) = U_hat_pathway_four(:, 1);   
    

    % % Commented on 12/08/2024 by Byeol
    % YS_Test_target_one_pathway_one = Ytest_target_one * V_pathway_one;
    % YS_Test_target_one_pathway_two = Ytest_target_one * V_pathway_two;
    % YS_Test_target_two_pathway_three = Ytest_target_two * V_pathway_three;
    % YS_Test_target_two_pathway_four = Ytest_target_two * V_pathway_four;
    % 
    % XS_Test_source_one_pathway_one = Xtest_source_one * Z_pathway_one;
    % XS_Test_source_two_pathway_two = Xtest_source_two * Z_pathway_two;
    % XS_Test_source_one_pathway_three = Xtest_source_one * Z_pathway_three;
    % XS_Test_source_two_pathway_four = Xtest_source_two * Z_pathway_four;
    % 
    % stats.latent_timeseries_source(indices==k, 1) = XS_Test_source_one_pathway_one(:, 1);
    % stats.latent_timeseries_source(indices==k, 2) = XS_Test_source_two_pathway_two(:, 1);
    % stats.latent_timeseries_source(indices==k, 3) = XS_Test_source_one_pathway_three(:, 1);
    % stats.latent_timeseries_source(indices==k, 4) = XS_Test_source_two_pathway_four(:, 1);
    % 
    % stats.latent_timeseries_target(indices==k, 1) = YS_Test_target_one_pathway_one(:, 1);
    % stats.latent_timeseries_target(indices==k, 2) = YS_Test_target_one_pathway_two(:, 1);
    % stats.latent_timeseries_target(indices==k, 3) = YS_Test_target_two_pathway_three(:, 1);
    % stats.latent_timeseries_target(indices==k, 4) = YS_Test_target_two_pathway_four(:, 1);


    % === Calculate latent correlation: Corr(T_hat, U_hat) ===

    % Within-pathways (On-target correlation): find correlation btw T and U
    % based on the optimized pair
    latent_correlation_pathway_one(k,:)=diag(corr(T_hat_pathway_one, U_hat_pathway_one));
    latent_correlation_pathway_four(k,:)=diag(corr(T_hat_pathway_four, U_hat_pathway_four));
    
    % Cross-pathways (Off-target correlation): find correlation btw T and U, 
    % but T was obtained from another pathway that has different source region.
    latent_correlation_pathway_one_crossed(k,:)=diag(corr(T_hat_pathway_two, U_hat_pathway_one)); % T_hat from X2-Y1 pathway & U_hat from X1-Y1 pathway
    latent_correlation_pathway_four_crossed(k,:)=diag(corr(T_hat_pathway_three, U_hat_pathway_four)); % T_hat from X1-Y2 pathway & U_hat from X2-Y2 pathway
    
end % End of Cross-validation loop 

% === Cross-validation results ===
% The latent correlations for within- and cross-pathways are stored in columns:
% Column 1: T_hat from X1-Y1 pathway and U_hat from X1-Y1 pathway [On, Within-pathway_one] 
% Column 2: T_hat from X2-Y1 pathway and U_hat from X1-Y1 pathway [Off, Cross-pathway_one]
% Column 3: T_hat from X1-Y2 pathway and U_hat from X2-Y2 pathway [Off, Cross-pathway_four]
% Column 4: T_hat from X1-Y1 pathway and U_hat from X2-Y2 pathway [On, Within-pathway_four]
stats.latent_correlations = [latent_correlation_pathway_one(:,1) latent_correlation_pathway_one_crossed(:,1) latent_correlation_pathway_four_crossed(:,1) latent_correlation_pathway_four(:,1) ];


% === Statistical tests ===
% Conduct t-tests on latent correlation coefficients for specific contrasts.
% Test whether on-target pathways are more functionally connected than off-target pathways.

% Interaction contrast [1 -1 -1 1], on-target > off-target
[~,p,~,stats.simple_correlation_interaction_ttest]=ttest(atanh(stats.simple_correlations(:,1))-atanh(stats.simple_correlations(:,2))-atanh(stats.simple_correlations(:,3))+atanh(stats.simple_correlations(:,4)));
stats.simple_correlation_interaction_ttest.p=p; % simple region average analysis
[~,p,~,stats.latent_correlation_interaction_ttest]=ttest(atanh(stats.latent_correlations(:,1))-atanh(stats.latent_correlations(:,2))-atanh(stats.latent_correlations(:,3))+atanh(stats.latent_correlations(:,4)));
stats.latent_correlation_interaction_ttest.p=p; % MPathI

% Pathway-specific contrasts: Pathway one (X1-Y1)
% Pathway one: [1 -1 0 0] on-target vs. off-target
[~,p,~,stats.simple_correlation_pathway_one_ttest]=ttest(atanh(stats.simple_correlations(:,1))-atanh(stats.simple_correlations(:,2)));
stats.simple_correlation_pathway_one_ttest.p=p; % simple region average analysis
[~,p,~,stats.latent_correlation_pathway_one_ttest]=ttest(atanh(stats.latent_correlations(:,1))-atanh(stats.latent_correlations(:,2)));
stats.latent_correlation_pathway_one_ttest.p=p; % MPathI

% Pathway two: [0 0 -1 1] on-target vs. off-target
[~,p,~,stats.simple_correlation_pathway_two_ttest]=ttest(atanh(stats.simple_correlations(:,4))-atanh(stats.simple_correlations(:,3)));
stats.simple_correlation_pathway_two_ttest.p=p; % simple region average analysis
[~,p,~,stats.latent_correlation_pathway_two_ttest]=ttest(atanh(stats.latent_correlations(:,4))-atanh(stats.latent_correlations(:,3)));
stats.latent_correlation_pathway_two_ttest.p=p; % MPathI

% on average, are 'on target' pathways more functionally connected than 'off target' pathways

%% Model parameters based on whole sample

% Perform hyperalignment across subjects

if do_alignment
    
    % Convert data into cell arrays, with one cell per subject
    % Each cell contains a Voxels x Time (Images) matrix
    for k=1:max(indices)
        to_align_dat_one{k}=source_one_obj.dat(:,indices==k);
        to_align_dat_two{k}=source_two_obj.dat(:,indices==k);
        to_align_dat_three{k}=target_one_obj.dat(:,indices==k);
        to_align_dat_four{k}=target_two_obj.dat(:,indices==k);
    end
    
    % Call Manning hyperalign.m
    % (Currently hyperaligns rows)
    aligned_dat_one = hyperalign(to_align_dat_one{:});
    aligned_dat_two = hyperalign(to_align_dat_two{:});
    aligned_dat_three = hyperalign(to_align_dat_three{:});
    aligned_dat_four = hyperalign(to_align_dat_four{:});
    
    for k=1:max(indices)
        source_one_obj.dat(:,indices==k)=aligned_dat_one{k};
        source_two_obj.dat(:,indices==k)=aligned_dat_two{k};
        target_one_obj.dat(:,indices==k)=aligned_dat_three{k};
        target_two_obj.dat(:,indices==k)=aligned_dat_four{k};
    end
end

% Fit an overall model for each pathway using the complete dataset, including both training and test data.
% xs: T, latent time series of the source region (X) [time x 1].
% ys: U, latent time series of the target region (Y) [time x 1].
% xl: loading/coefficient for the source region (X), calculated for sign checking.
[xl_pathway_one, ~, xs_pathway_one, ys_pathway_one] = plsregress(source_one_obj.dat',target_one_obj.dat', ndim);
[xl_pathway_two, ~, xs_pathway_two, ys_pathway_two] = plsregress(source_two_obj.dat',target_one_obj.dat', ndim);
[xl_pathway_three, ~, xs_pathway_three, ys_pathway_three] = plsregress(source_one_obj.dat',(target_two_obj.dat)', ndim);
[xl_pathway_four, ~, xs_pathway_four, ys_pathway_four] = plsregress(source_two_obj.dat',(target_two_obj.dat)', ndim);

% Optionally flip the latent variable signs to ensure consistency
if flip_maps
    % Determine the sign flip based on spatial correlation between latent X loadings and source region averages
    sign_pathway(1) = sign(corr(mean(xl_pathway_one, 2), mean(source_one_obj.dat, 2)));
    sign_pathway(2) = sign(corr(mean(xl_pathway_two, 2), mean(source_two_obj.dat, 2)));
    sign_pathway(3) = sign(corr(mean(xl_pathway_three, 2), mean(source_one_obj.dat, 2)));    
    sign_pathway(4) = sign(corr(mean(xl_pathway_four,2), mean(source_two_obj.dat,2)));
    
    % Adjust scores for X (source) and Y (target) latent variables based on sign flips
    ys_pathway_one = ys_pathway_one * sign_pathway(1);
    ys_pathway_two = ys_pathway_two * sign_pathway(2);
    ys_pathway_three = ys_pathway_three * sign_pathway(3);
    ys_pathway_four = ys_pathway_four * sign_pathway(4);
    
    xs_pathway_one = xs_pathway_one * sign_pathway(1);
    xs_pathway_two = xs_pathway_two * sign_pathway(2);
    xs_pathway_three = xs_pathway_three * sign_pathway(3);
    xs_pathway_four = xs_pathway_four * sign_pathway(4);
  
    % Adjust overall latent time series to match the sign flip
   for p = 1:4
        stats.latent_timeseries_source(:, p) = stats.latent_timeseries_source(:, p) * sign_pathway(p);
        stats.latent_timeseries_target(:, p) = stats.latent_timeseries_target(:, p) * sign_pathway(p);
    end
    stats.pathway_sign = sign_pathway;
    
%     stats.latent_timeseries(:, 1) = stats.latent_timeseries(:, 1) * sign(corr(mean(xl_pathway_one,2),mean(source_one_obj.dat,2)));
%     stats.latent_timeseries(:, 2) = stats.latent_timeseries(:, 2) * sign(corr(mean(xl_pathway_four,2),mean(source_two_obj.dat,2)));

end

% added by Tor, 2/23/2025
stats.path1_overall_xval_r = corr(stats.latent_timeseries_source(:, 1), stats.latent_timeseries_target(:, 1));
stats.path2_overall_xval_r = corr(stats.latent_timeseries_source(:, 4), stats.latent_timeseries_target(:, 4));

stats.path1_overall_xval_dot = dot(stats.latent_timeseries_source(:, 1), stats.latent_timeseries_target(:, 1));
stats.path2_overall_xval_dot = dot(stats.latent_timeseries_source(:, 4), stats.latent_timeseries_target(:, 4));

% added by Byeol on 12/12/2024
stats.T_pathway_one = xs_pathway_one;
stats.U_pathway_one = ys_pathway_one;

stats.T_pathway_four = xs_pathway_four;
stats.U_pathway_four = ys_pathway_four;

% Compute latent weights (Z) for source regions and pattern weights (V) for
% target regions with updated T and U.
% Reminder: Z = pinv(X) * U, V = pinv(Y) * T [Voxels x 1]
stats.Z_pathway_one = pinv(source_one_obj.dat') * ys_pathway_one; 
stats.V_pathway_one = pinv(target_one_obj.dat') * xs_pathway_one; 

stats.Z_pathway_four = pinv(source_two_obj.dat') * ys_pathway_four;
stats.V_pathway_four = pinv(target_two_obj.dat') * xs_pathway_four;


% Perform bootstrap resampling if specified
if do_boot    
    bs_source_one_dat=source_one_obj.dat';
    bs_source_two_dat=source_two_obj.dat';
    bs_target_one_dat=target_one_obj.dat';
    bs_target_two_dat=target_two_obj.dat';
    for i=1:nboot
        % Resample indices with replacement at the fold level
        [~,rand_subs]=datasample(1:max(indices), max(indices)); %randomly replace whole blocks with replacement
        count_vec=1:length(indices);

        rand_inds=[];
        for ii=1:max(indices)
            rand_inds=[rand_inds count_vec(indices==rand_subs(ii))];
        end

        % Compute bootstrap estimates for source and target weights
        bs_V_pathway_one(i,:)= bootPLS_target_pattern_weights(bs_source_one_dat(rand_inds,:),bs_target_one_dat(rand_inds,:),ndim,flip_maps);
        bs_Z_pathway_one(i,:)= bootPLS_source_pattern_weights(bs_source_one_dat(rand_inds,:),bs_target_one_dat(rand_inds,:),ndim,flip_maps);
        
        bs_V_pathway_four(i,:)= bootPLS_target_pattern_weights(bs_source_two_dat,bs_target_two_dat,ndim,flip_maps);
        bs_Z_pathway_four(i,:)= bootPLS_source_pattern_weights(bs_source_two_dat,bs_target_two_dat,ndim,flip_maps);
    end
    
    
%     bs_V_pathway_one= bootstrp(nboot,@bootPLS_target_pattern_weights,source_one_obj.dat',target_one_obj.dat',ndim,flip_maps);
%     bs_Z_pathway_one= bootstrp(nboot,@bootPLS_source_pattern_weights,source_one_obj.dat',target_one_obj.dat',ndim,flip_maps);
   
%     bs_V_pathway_four= bootstrp(nboot,@bootPLS_target_pattern_weights,source_two_obj.dat',target_two_obj.dat',ndim,flip_maps);
%     bs_Z_pathway_four= bootstrp(nboot,@bootPLS_source_pattern_weights,source_two_obj.dat',target_two_obj.dat',ndim,flip_maps);
    
    % Compute bootstrap statistics for Z and V (weights)
    bs_Z=(mean(bs_Z_pathway_one)./std(bs_Z_pathway_one))';
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    bs_stat=statistic_image;
    bs_stat.volInfo=source_one_obj.volInfo;
    bs_stat.dat=stats.Z_pathway_one;
    bs_stat.p=bs_P;
    bs_stat.removed_voxels=source_one_obj.removed_voxels;
    stats.PLS_bootstrap_stats_Z(1)=bs_stat;
    
    bs_Z=(mean(bs_Z_pathway_four)./std(bs_Z_pathway_four))';
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    bs_stat=statistic_image;
    bs_stat.volInfo=source_two_obj.volInfo;
    bs_stat.dat=stats.Z_pathway_four;
    bs_stat.p=bs_P;
    bs_stat.removed_voxels=source_two_obj.removed_voxels;
    stats.PLS_bootstrap_stats_Z(2)=bs_stat;
    
    
    bs_Z=(mean(bs_V_pathway_one)./std(bs_V_pathway_one))';
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    bs_stat=statistic_image;
    bs_stat.volInfo=target_one_obj.volInfo;
    bs_stat.dat=stats.V_pathway_one;
    bs_stat.p=bs_P;
    bs_stat.removed_voxels=target_one_obj.removed_voxels;
    stats.PLS_bootstrap_stats_V(1)=bs_stat;
    
    bs_Z=(mean(bs_V_pathway_four)./std(bs_V_pathway_four))';
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    bs_stat=statistic_image;
    bs_stat.volInfo=target_two_obj.volInfo;
    bs_stat.dat=stats.V_pathway_four;
    bs_stat.p=bs_P;
    bs_stat.removed_voxels=target_two_obj.removed_voxels;
    stats.PLS_bootstrap_stats_V(2)=bs_stat;
end


%% Output data objects
if do_roi

    stats.source_one_obj=source_one_obj;
    stats.source_two_obj=source_two_obj;
    stats.target_one_obj=target_one_obj;
    stats.target_two_obj=target_two_obj;

end
%% Add results report with narrative text
stats = add_results_report(stats);

%% Optional plot
if do_plot

    colors = {[.3 .4 1] [.4 .5 .6] [.6 .5 .4] [1 .7 .3]};

    % Plot the correlations among region averages:
    disp('Simple ROI correlations')
    create_figure('Simple ROI correlations');

    barplot_columns(stats.simple_correlations, 'names', pathwaynames, 'colors', colors, 'nofigure');
    title('Simple ROI correlations')

    % Plot the correlations among PLS-optimized latent timeseries:
    disp('MPathI correlations')
    create_figure('MPathI correlations');

    barplot_columns(stats.latent_correlations, 'names', pathwaynames, 'colors', colors, 'nofigure');
    title('MPathI correlations')

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


function stats = add_results_report(stats)


% Latent MPathI correlations

r = stats.latent_correlations;
[h, p, ci, tstat] = ttest(atanh(r));  % t-test on Fisher's z-transformed latent correlations
if any(isnan(h))
    warning('The latent correlations have values of 1 or -1, which do not seem correct. The following report is calculated after removing the rows (folds) where r = 1 or -1.');
    r(any(abs(r') == 1),:) = [];
    [h, p, ci, tstat] = ttest(atanh(r));  % t-test on Fisher's z-transformed latent correlations. Same as fisherz.m
end
rm = mean(r);
rse = ste(r);


sigstr = {'not' ''};

stats_report{1} = sprintf('Latent MPathI correlations were %3.2f (sd = %3.2f) for Pathway 1 and %3.2f (sd = %3.2f) for Pathway 2.', rm(1), rse(1), rm(4), rse(4));
stats_report{2} = sprintf('Correlations were %s significant for Pathway 1, t(%3.1f) = %3.2f, p = %3.6f, and %s significant for Pathway 2, t(%3.1f) = %3.2f, p = %3.6f.', ...
    sigstr{h(1)+1}, tstat.df(1), tstat.tstat(1), p(1), sigstr{h(4)+1}, tstat.df(4), tstat.tstat(4), p(4));

istat = stats.latent_correlation_interaction_ttest;

stats_report{3} = sprintf('The interaction between on-targed (expected) and off-target pathways was %s signficant, t(%3.1f) = %3.2f, p = %3.6f.', ...
    sigstr{(istat.p < 0.05) + 1}, istat.df(1), istat.tstat(1), istat.p(1));

istat = stats.latent_correlation_pathway_one_ttest;

stats_report{4} = sprintf('Pathway 1 on-target vs. off-target correlations were %s signficant, t(%3.1f) = %3.2f, p = %3.6f.', ...
    sigstr{(istat.p < 0.05) + 1}, istat.df(1), istat.tstat(1), istat.p(1));

istat = stats.latent_correlation_pathway_two_ttest;

stats_report{5} = sprintf('Pathway 2 on-target vs. off-target correlations were %s signficant, t(%3.1f) = %3.2f, p = %3.6f.', ...
    sigstr{(istat.p < 0.05) + 1}, istat.df(1), istat.tstat(1), istat.p(1));

stats_report{6} = '       ';

wh_last = length(stats_report);

% Simple correlations

r = stats.simple_correlations;
[h, p, ci, tstat] = ttest(atanh(r));  % t-test on Fisher's z-transformed latent correlations. Same as fisherz.m
if any(isnan(h))
    warning('The simple correlations have values of 1 or -1, which do not seem correct. The following report is calculated after removing the rows (folds) where r = 1 or -1.');
    r(any(abs(r') == 1),:) = [];
    [h, p, ci, tstat] = ttest(atanh(r));  % t-test on Fisher's z-transformed latent correlations. Same as fisherz.m
end
rm = mean(r);
rse = ste(r);

sigstr = {'not' ''};

stats_report{wh_last + 1} = sprintf('Simple ROI correlations were %3.2f (sd = %3.2f) for Pathway 1 and %3.2f (sd = %3.2f) for Pathway 2.', rm(1), rse(1), rm(4), rse(4));
stats_report{wh_last + 2} = sprintf('Correlations were %s significant for Pathway 1, t(%3.1f) = %3.2f, p = %3.6f, and %s significant for Pathway 2, t(%3.1f) = %3.2f, p = %3.6f.', ...
    sigstr{h(1)+1}, tstat.df(1), tstat.tstat(1), p(1), sigstr{h(4)+1}, tstat.df(4), tstat.tstat(4), p(4));

istat = stats.simple_correlation_interaction_ttest;

stats_report{wh_last + 3} = sprintf('The interaction between on-targed (expected) and off-target pathways was %s signficant, t(%3.1f) = %3.2f, p = %3.6f.', ...
    sigstr{(istat.p < 0.05) + 1}, istat.df(1), istat.tstat(1), istat.p(1));

istat = stats.simple_correlation_pathway_one_ttest;

stats_report{wh_last + 4} = sprintf('Pathway 1 on-target vs. off-target correlations were %s signficant, t(%3.1f) = %3.2f, p = %3.6f.', ...
    sigstr{(istat.p < 0.05) + 1}, istat.df(1), istat.tstat(1), istat.p(1));

istat = stats.simple_correlation_pathway_two_ttest;

stats_report{wh_last + 5} = sprintf('Pathway 2 on-target vs. off-target correlations were %s signficant, t(%3.1f) = %3.2f, p = %3.6f.', ...
    sigstr{(istat.p < 0.05) + 1}, istat.df(1), istat.tstat(1), istat.p(1));

stats_report{wh_last + 6} = '       ';

stats.stats_report = char(stats_report{:});

end
