function stats = model_brain_pathway(obj,source_one,source_two,target_one,target_two,varargin)
% Models connections between a pair of brain regions using partial least squares regression
%
% Usage:
% ::
%
%    stats = model_brain_pathway(obj,source_one,source_two,target_one,target_two,varargin);
%
% - Default: 10-fold cross-validation (no grouping; see optional inputs) for testing connection strength
% - Default: block bootstrap with random selection of observations in blocks for significance of pattern weights
% - If running across multiple participants, enter 'Indices' (see below) to properly hangle subject grouping
%
% This is a method for an fmri_data object that attempts to identify a set of
% weights in one brain region (X, the source) that can be used to explain activity
% in another brain region (Y, the target region). This many to many mapping is
% performed using Partial Least Squares regression, and is based on the
% idea that connected neural populations in the two brain regions produce
% correlated signals (e.g., fMRI activation). This function performs the following steps:
%
% 1. Perform a typical cross-validated connectivity analysis as a comparison condition
%   - estimate simple correlations using k-fold CV
%   - 4 paths connecting 2 regions, using region averages:
%   - X1->Y1, X2->Y1, X1->Y2, X2->Y2
%   - These are described as pathways 1-4.
%   - Custom holdout set is specified by a vector of counts, 1 per holdout set
%
% 2. Cross-validate PLS model estimating pathways
%   - if specified, use input cross-validation indices
%   - if specified, create hyperalignment model on training data, and warp
%   test data into group mean space
%   - estimate PLS models for the four pathways described above, with one
%   latent source
%   - for all pathways, estimate V weights by regressing voxels in Y
%   onto latent timeseries for X and estimate Z weights by regressing
%   voxels in X onto latent timeseries for Y (estimated on training data)
%   - apply Z and V weights to test data to get estimate of latent
%   population activity (i.e., latent timeseries) and compute the 
%   correlation between timeseries using Pearson's correlation coefficient.
%   This is done between Yhat_V1 and Y1, Yhat_V2 and Y1, Yhat_V1 and Y2, 
%   and Yhat_V2 and Y2, with expectation that crossing the pathways leads 
%   to lower correlations between Yhat_V2 and Y1 and Yhat_V1 and Y2.
%
% 3. Estimate the PLS models on the full sample, 
%    - hyperalign if requested
%    - bootstrap for voxel-level inference on V and Z weights if requested
%
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
%
% :Inputs:
%
%   **obj:**
%        An fmri_data image object with one or more images loaded (e.g., single-trial
%        or preprocessed time-series)
%
%   **source_one:**
%        An fmri_data object whose .dat field is a binary mask specifying
%        voxels that belong to source 1 (X in PLS), or X1
%
%   **source_two:**
%        An fmri_data object whose .dat field is a binary mask specifying
%        voxels that belong to source 2 (X in PLS), or X2
%
%   **target_one:**
%        An fmri_data object whose .dat field is a binary mask specifying
%        voxels that belong to target 1 (Y in PLS), or Y1
%
%   **target_two:**
%        An fmri_data object whose .dat field is a binary mask specifying
%        voxels that belong to target 2 (Y in PLS), or Y2
%
%
% :Optional inputs:
%
%   **'nboot'**
%      Followed by the number of bootstrap samples to conduct for
%      estimating voxel-wise significance of V and Z weights
%
% **'Indices'**
%      Followed by an integer vector specifying participant (or blocking) ID
%      Used in two ways:
%      (1) to define cross-validation holdout set for each observation, n_obs x 1. This is useful if you want to,
%      e.g., leave out all images from a participant.
%      (2) to define blocks for block bootstrap
%
%   **'Align'**
%      Perform hyperalignment, requires input 'Indices' to be specified
%      with each index corresponding to a subject
%
% :Outputs:
%
%   **stats:**
%        Structure including:
%           - simple_correlations:  correlations among region averages,
%             order: 1->1, 1->2, 2->1, 2->2
%           - latent_correlations:  correlations among PLS-optimized latent timeseries,
%             order: 1->1, 1->2, 2->1, 2->2
%           - latent_timeseries:
%             latent_correlation_interaction_ttest: T-test on Fisher-transformed on-target vs. off-target
%             simple_correlation_interaction_ttest: [1×1 struct]
%             latent_correlation_pathway_one_ttest: T-test on
%             Fisher-transformed, pathway 1. Pathway 1 connection is
%             non-zero.
%             latent_correlation_pathway_two_ttest: T-test on
%             Fisher-transformed, pathway 2. Pathway 2 connection is
%             non-zero.
%             simple_correlation_pathway_one_ttest: [1×1 struct]
%             simple_correlation_pathway_two_ttest: [1×1 struct]
%           - PLS_betas (not included now): fmri_data object with PLS regression coefficients
%           - PLS_bootstrap_stats_Z: statistic image for PLS regression coefficients
%           - PLS_bootstrap_stats_V: statistic image for PLS regression coefficients
%           - Z_pathway_one/four: estimated patterns in source regions that covary with each target 
%           - V_pathway_one/four: estimated patterns in target regions that covary with each source
%           - source_one_obj: resampled and masked fmri data object for source one
%           - source_two_obj: resampled and masked fmri data object for source two
%           - target_one_obj: resampled and masked fmri data object for target one
%           - target_two_obj: resampled and masked fmri data object for target one
%
%
% :Examples:
% ::
%
% % Load multi-subject image object (toy example, 6 images per subject)
% imgs = load_image_set('bmrk3'); 
%
% % Define ROIs
% vpl = select_atlas_subset(atlas_obj, {'VPL'}, 'flatten');
% pbn = select_atlas_subset(atlas_obj, {'pbn'}, 'flatten');
% orthviews(vpl);
% orthviews(pbn);
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
% figure; barplot_columns(stats.simple_correlations, 'names', {'vpl->dpins' 'vpl->cea' 'pbn->dpins' 'pbn->cea'});
%
% % Plot the correlations among PLS-optimized latent timeseries:
% figure; barplot_columns(stats.latent_correlations, 'names', {'vpl->dpins' 'vpl->cea' 'pbn->dpins' 'pbn->cea'});
%
% % Run it again with bootstrapping too:
% stats = model_brain_pathway(imgs, vpl, pbn, dpins, cea, 'Indices', indx, 'nboot', 1000);
%
% :See also:
%

% ..
%    Programmers' notes:
% 10/10/2019 Wrote funcion (Phil Kragel)
% 4/14/2020 Major overhaul and documentation (Phil Kragel and Tor Wager)
% 4/15/2020 Add more documentation (Phil Kragel)
% 9/16/2020 Tor Wager - update documentation
%
% Notes for future development:
% - could use xval_stratified_holdout_leave_whole_subject_out in future versions for cross-validation
% - could write separate hyperalignment object method
% 
%
% ..

%% Get defaults and initialize user inputs
if any(strcmp(varargin,'Align'))
    do_alignment=true;
else
    do_alignment=false;
end

if any(strcmp(varargin,'nboot'))
    do_boot=true;
    nboot  = varargin{find(strcmp(varargin,'nboot'))+1};
    
else
    do_boot=false;
end

indices=crossvalind('Kfold',size(obj.dat,2),10); %do 10-fold by default.
if any(strcmp(varargin,'Indices'))
    indices = varargin{find(strcmp(varargin,'Indices'))+1};
end

flip_maps=true;


ndim=1; %number of latent dims for PLS

%% convert data and masks to masked objects

source_one_obj=apply_mask(obj,source_one);
source_two_obj=apply_mask(obj,source_two);
target_one_obj=apply_mask(obj,target_one);
target_two_obj=apply_mask(obj,target_two);


%% Perform a typical cross-validated connectivity analysis as a comparison condition
% estimate simple correlations using k-fold CV
% 4 paths connecting 2 regions, using region averages:
% X1<->Y1, X2<->Y1, X1<->Y2, X2<->Y2
% These are pathways 1-4.
% Custom holdout set is specified by a vector of integers, 1 per holdout set

for k=1:max(indices) %for each fold
    
    train = indices~=k; %train on all but that fold
    test = ~train; %test on that fold
    
    %simple linear association of means for comparison
    beta_pathway_one_mean = glmfit(mean(source_one_obj.dat(:,train))',mean(target_one_obj.dat(:,train))');
    beta_pathway_two_mean = glmfit(mean(source_two_obj.dat(:,train))',mean(target_one_obj.dat(:,train))');
    beta_pathway_three_mean = glmfit(mean(source_one_obj.dat(:,train))',mean(target_two_obj.dat(:,train))');
    beta_pathway_four_mean = glmfit(mean(source_two_obj.dat(:,train))',mean(target_two_obj.dat(:,train))');
    
    
    %simple regression model predictions
    yhat_pathway_one_mean(test,:)=[ones(size(mean(source_one_obj.dat(:,test))',1),1) mean(source_one_obj.dat(:,test))'] * beta_pathway_one_mean;
    yhat_pathway_two_mean(test,:)=[ones(size(mean(source_two_obj.dat(:,test))',1),1) mean(source_two_obj.dat(:,test))'] * beta_pathway_two_mean;
    yhat_pathway_three_mean(test,:)=[ones(size(mean(source_one_obj.dat(:,test))',1),1) mean(source_one_obj.dat(:,test))'] * beta_pathway_three_mean;
    yhat_pathway_four_mean(test,:)=[ones(size(mean(source_two_obj.dat(:,test))',1),1) mean(source_two_obj.dat(:,test))'] * beta_pathway_four_mean;
    
    %simple correlations of means
    
    
    stats.simple_correlations(k,:)=[corr(yhat_pathway_one_mean(test),mean(target_one_obj.dat(:,test))') corr(yhat_pathway_one_mean(test),mean(target_two_obj.dat(:,test))') corr(yhat_pathway_four_mean(test),mean(target_one_obj.dat(:,test))') corr(yhat_pathway_four_mean(test),mean(target_two_obj.dat(:,test))')];
    
end


%% estimate pattern-based connectivity with kfold or user-specified xval

for k=1:max(indices)
    
    un_inds=unique(indices);
  
    
    % Convert to cell array, one cell per subject
    % Voxels x Images matrix
    
    for kk=1:max(indices)
        to_align_dat_one{kk}=source_one_obj.dat(:,indices==kk);
        to_align_dat_two{kk}=source_two_obj.dat(:,indices==kk);
        to_align_dat_three{kk}=target_one_obj.dat(:,indices==kk);
        to_align_dat_four{kk}=target_two_obj.dat(:,indices==kk);
        
    end
    
    if do_alignment
        % Call Manning hyperalign.m
        % (Currently hyperaligns rows)
             
        aligned_dat_one = hyperalign(to_align_dat_one{un_inds~=k});
        aligned_dat_two = hyperalign(to_align_dat_two{un_inds~=k});
        aligned_dat_three = hyperalign(to_align_dat_three{un_inds~=k});
        aligned_dat_four = hyperalign(to_align_dat_four{un_inds~=k});
        
        source_one_train=[aligned_dat_one{:}];
        source_two_train=[aligned_dat_two{:}];
        target_one_train=[aligned_dat_three{:}];
        target_two_train=[aligned_dat_four{:}];
        
        [~, source_one_test] = procrustes(mean(cat(3,to_align_dat_one{un_inds~=k}),3), to_align_dat_one{un_inds==k});
        [~, source_two_test] = procrustes(mean(cat(3,to_align_dat_two{un_inds~=k}),3), to_align_dat_two{un_inds==k});
        [~, target_one_test] = procrustes(mean(cat(3,to_align_dat_three{un_inds~=k}),3), to_align_dat_three{un_inds==k});
        [~, target_two_test] = procrustes(mean(cat(3,to_align_dat_four{un_inds~=k}),3), to_align_dat_four{un_inds==k});
        
    else
        
        source_one_train=[to_align_dat_one{un_inds~=k}];
        source_two_train=[to_align_dat_two{un_inds~=k}];
        target_one_train=[to_align_dat_three{un_inds~=k}];
        target_two_train=[to_align_dat_four{un_inds~=k}];
        
        source_one_test=[to_align_dat_one{un_inds==k}];
        source_two_test=[to_align_dat_two{un_inds==k}];
        target_one_test=[to_align_dat_three{un_inds==k}];
        target_two_test=[to_align_dat_four{un_inds==k}];
    end
    
    
    
    
    %pls models to estimate correlation between latent sources
    [~,~,xs_pathway_one,ys_pathway_one] = plsregress(source_one_train',target_one_train',ndim);
    [~,~,xs_pathway_two,ys_pathway_two] = plsregress(source_two_train',target_one_train',ndim);
    [~,~,xs_pathway_three,ys_pathway_three] = plsregress(source_one_train',target_two_train',ndim);
    [~,~,xs_pathway_four,ys_pathway_four] = plsregress(source_two_train',target_two_train',ndim);
    
 
     %
    Z_pathway_one = pinv(source_one_train') * ys_pathway_one;
    V_pathway_one = pinv(target_one_train') * xs_pathway_one;
    
    Z_pathway_two = pinv(source_two_train') * ys_pathway_two;
    V_pathway_two = pinv(target_one_train') * xs_pathway_two;
    
    Z_pathway_three = pinv(source_one_train') * ys_pathway_three;
    V_pathway_three = pinv(target_two_train') * xs_pathway_three;
    
    Z_pathway_four = pinv(source_two_train') * ys_pathway_four;
    V_pathway_four = pinv(target_two_train') * xs_pathway_four;
    
    
    
    %correlation of latent dimensions
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
    
    %optimized models
    
    YS_Test_target_one_pathway_one = Ytest_target_one*V_pathway_one;
    YS_Test_target_one_pathway_two = Ytest_target_one*V_pathway_two;
    YS_Test_target_two_pathway_three = Ytest_target_two*V_pathway_three;
    YS_Test_target_two_pathway_four = Ytest_target_two*V_pathway_four;
    
    
    XS_Test_source_one_pathway_one = Xtest_source_one*Z_pathway_one;
    XS_Test_source_two_pathway_two = Xtest_source_two*Z_pathway_two;
    XS_Test_source_one_pathway_three = Xtest_source_one*Z_pathway_three;
    XS_Test_source_two_pathway_four = Xtest_source_two*Z_pathway_four;
    
    stats.latent_timeseries(indices==k,1)=YS_Test_target_one_pathway_one(:,1);
    stats.latent_timeseries(indices==k,2)=YS_Test_target_two_pathway_four(:,1);
    
    %optimized pathways
    latent_correlation_pathway_one(k,:)=diag(corr(YS_Test_target_one_pathway_one,XS_Test_source_one_pathway_one));
    latent_correlation_pathway_four(k,:)=diag(corr(YS_Test_target_two_pathway_four,XS_Test_source_two_pathway_four));
    
    %switched targets
    latent_correlation_pathway_one_crossed(k,:)=diag(corr(YS_Test_target_one_pathway_two,XS_Test_source_one_pathway_one));%;XS_Test_source_one_pathway_three
    latent_correlation_pathway_four_crossed(k,:)=diag(corr(YS_Test_target_two_pathway_three,XS_Test_source_two_pathway_four));%XS_Test_source_two_pathway_two
    
    
    
    
end


stats.latent_correlations=[latent_correlation_pathway_one(:,1) latent_correlation_pathway_one_crossed(:,1) latent_correlation_pathway_four_crossed(:,1) latent_correlation_pathway_four(:,1) ];

% on average, are 'on target' pathways more functionally connected than 'off target' pathways
[~,p,~,stats.latent_correlation_interaction_ttest]=ttest(atanh(stats.latent_correlations(:,1))-atanh(stats.latent_correlations(:,2))-atanh(stats.latent_correlations(:,3))+atanh(stats.latent_correlations(:,4)));
stats.latent_correlation_interaction_ttest.p=p;
[~,p,~,stats.simple_correlation_interaction_ttest]=ttest(atanh(stats.simple_correlations(:,1))-atanh(stats.simple_correlations(:,2))-atanh(stats.simple_correlations(:,3))+atanh(stats.simple_correlations(:,4)));
stats.simple_correlation_interaction_ttest.p=p;

[~,p,~,stats.latent_correlation_pathway_one_ttest]=ttest(atanh(stats.latent_correlations(:,1))-atanh(stats.latent_correlations(:,2)));
stats.latent_correlation_pathway_one_ttest.p=p;

[~,p,~,stats.latent_correlation_pathway_two_ttest]=ttest(atanh(stats.latent_correlations(:,4))-atanh(stats.latent_correlations(:,3)));
stats.latent_correlation_pathway_two_ttest.p=p;

[~,p,~,stats.simple_correlation_pathway_one_ttest]=ttest(atanh(stats.simple_correlations(:,1))-atanh(stats.simple_correlations(:,2)));
stats.simple_correlation_pathway_one_ttest.p=p;

[~,p,~,stats.simple_correlation_pathway_two_ttest]=ttest(atanh(stats.simple_correlations(:,4))-atanh(stats.simple_correlations(:,3)));
stats.simple_correlation_pathway_two_ttest.p=p;


%% Model parameters based on whole sample

% Perform hyperalignment


if do_alignment
    
    % Convert to cell array, one cell per subject
    % Voxels x Images matrix
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

[xl_pathway_one,~,xs_pathway_one,ys_pathway_one] = plsregress(source_one_obj.dat',target_one_obj.dat',ndim);
[xl_pathway_four,~,xs_pathway_four,ys_pathway_four] = plsregress(source_two_obj.dat',(target_two_obj.dat)',ndim);

if flip_maps
ys_pathway_one=ys_pathway_one*sign(corr(mean(xl_pathway_one,2),mean(source_one_obj.dat,2)));
ys_pathway_four=ys_pathway_four*sign(corr(mean(xl_pathway_four,2),mean(source_two_obj.dat,2)));

xs_pathway_one=xs_pathway_one*sign(corr(mean(xl_pathway_one,2),mean(source_one_obj.dat,2)));
xs_pathway_four=xs_pathway_four*sign(corr(mean(xl_pathway_four,2),mean(source_two_obj.dat,2)));

stats.latent_timeseries(:,1)=stats.latent_timeseries(:,1)*sign(corr(mean(xl_pathway_one,2),mean(source_one_obj.dat,2)));
stats.latent_timeseries(:,2)=stats.latent_timeseries(:,2)*sign(corr(mean(xl_pathway_four,2),mean(source_two_obj.dat,2)));
end

stats.Z_pathway_one = pinv(source_one_obj.dat') * ys_pathway_one;
stats.V_pathway_one = pinv(target_one_obj.dat') * xs_pathway_one;

stats.Z_pathway_four = pinv(source_two_obj.dat') * ys_pathway_four;
stats.V_pathway_four = pinv(target_two_obj.dat') * xs_pathway_four;

if do_boot
    
    bs_source_one_dat=source_one_obj.dat';
    bs_source_two_dat=source_two_obj.dat';
    bs_target_one_dat=target_one_obj.dat';
    bs_target_two_dat=target_two_obj.dat';
    for i=1:nboot
        
        [~,rand_subs]=datasample(1:max(indices), max(indices)); %randomly replace whole blocks with replacement
        count_vec=1:length(indices);
        for ii=1:max(indices)
            
            rand_inds(indices==ii)=count_vec(indices==rand_subs(ii));
            
        end
        
        bs_V_pathway_one(i,:)= bootPLS_target_pattern_weights(bs_source_one_dat(rand_inds,:),bs_target_one_dat(rand_inds,:),ndim,flip_maps);
        bs_Z_pathway_one(i,:)= bootPLS_source_pattern_weights(bs_source_one_dat(rand_inds,:),bs_target_one_dat(rand_inds,:),ndim,flip_maps);
        
        
        bs_V_pathway_four(i,:)= bootPLS_target_pattern_weights(bs_source_two_dat,bs_target_two_dat,ndim,flip_maps);
        bs_Z_pathway_four(i,:)= bootPLS_source_pattern_weights(bs_source_two_dat,bs_target_two_dat,ndim,flip_maps);
    end
    
    
%     bs_V_pathway_one= bootstrp(nboot,@bootPLS_target_pattern_weights,source_one_obj.dat',target_one_obj.dat',ndim,flip_maps);
%     bs_Z_pathway_one= bootstrp(nboot,@bootPLS_source_pattern_weights,source_one_obj.dat',target_one_obj.dat',ndim,flip_maps);
    
    
%     bs_V_pathway_four= bootstrp(nboot,@bootPLS_target_pattern_weights,source_two_obj.dat',target_two_obj.dat',ndim,flip_maps);
%     bs_Z_pathway_four= bootstrp(nboot,@bootPLS_source_pattern_weights,source_two_obj.dat',target_two_obj.dat',ndim,flip_maps);
    
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


%% output data objects
stats.source_one_obj=source_one_obj;
stats.source_two_obj=source_two_obj;
stats.target_one_obj=target_one_obj;
stats.target_two_obj=target_two_obj;

end


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