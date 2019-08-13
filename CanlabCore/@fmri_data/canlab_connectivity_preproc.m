function [dat, roi_val, maskdat, beta_dat, beta_roi_val] = canlab_connectivity_preproc(dat, varargin)

% This function prepares data for connectivity analysis by removing nuisance
% variables and temporal filtering (high, low, or bandpass filter). This also
% can extract values from given masks and return averaged activity or pattern
% expression values.
%
% :Usage:
% ::
%
%    [preprocessed_dat, roi_val, maskdat] = canlab_connectivity_preproc(dat, varargin)
%    [preprocessed_dat, roi_val, maskdat, (optional) beta_dat, beta_roi_val] = canlab_connectivity_preproc(dat, varargin)
%
% :Features:
%
%   - can regress out nuisance variables with any additional nuisance matrix
%   - can remove signal from ventricle and white matter (default: use
%     canonical ventricle and white matter masks)
%   - can do temporal filtering, including high-pass, low-pass, or bandpass
%     filtering (it uses conn_filter.m from conn toolbox, which is using
%     fast fourier transform; see subfunction below)
%   - can extract data from given ROIs, and return averaged value or pattern
%     expression value (dot-product).
%   - Subsequent nuisance regression step is orthogonalized with respect to
%     pass-filter, consistent with principles in Lindquist 2019
%   - It is possible to extract ROI averages or linear patterns, save them,
%     and run this function afterwards. This will give the same result as
%     filtering each  voxel and then calculating the ROI average/pattern
%     response.
%   - (optional) can run additional GLM model with additional regressors if they are 
%     specified. This runs the GLM *additionally*, therefore do not change
%     residualized image data or roi_vals. This GLM uses the same nuisance
%     matrix from above. output: beta_dat  
%
% *Steps in order [with defaults]:*
%   1. Remove nuisance covariates (and linear trend if requested)
%   2. Remove ventricle and white matter - default: uses canonical images
%   3. High/low/bandpass filter on BOTH data and covs before regression (plot)
%   4. Residualization using regression (plot)
%   5. (optional) Run additional GLM using the same additional preprocessing
%   6. Windsorize based on distribution of full data matrix (plot)
%   7. Extract region-by-region average ROI or pattern expression data
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2014  Wani Woo
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
%   **dat:**
%        fmri_data object with data
%
%   **dat.covariate:**
%        basic nuisance matrix
%
% :Optional Inputs:
%
%   **additional_nuisance**
%        When you have additional nuisance variables that you want
%        regress out from the data, you can use this option. This
%        option should be followed by a nuisance matrix (or values).
%        The matrix should have the same number of rows with the
%        number of images.
%
%   **vw**
%        When you want to regress out signals from ventricle and
%        white matter, you can use this option. Default is using canonical
%        white matter and ventricle images (but eroded versions). 
%
%        You can also choose what to use to remove ventricle and
%        white matter signal between raw data or top 5 PCA
%        components (default). You can just put 'raw' if you want to
%        use raw signal than PCA compoenents.
%        also see: canlab_extract_ventricle_wm_timeseries.m
%        canlab_create_wm_ventricle_masks.m)
%        - *Example:* 'vw','raw'
%
%   **windsorize:**
%        Windsorizing entire data matrix to k x STD.
%        - *Example:* 'windsorize', 5 (windsorize to 5 STD)
%
%   **linear_trend:**
%        This option will include the linear trend to nuisance variables.
%
%   **hpf', 'lpf', or 'bpf:**
%        This option will do temporal filtering.
%        - 'hpf': high pass filter. This option should be followed by
%                 the lower bound of the frequency (e.g., .01 Hz [= 100
%                 sec]), then TR
%        - 'lpf': low pass filter. This option should be followed by
%                 the upper bound of the frequency (e.g., .25 Hz [= 4 sec]), then TR.
%        - 'bpf': bandpass filter. This should be followed by lower
%                 and upper bounds of the frequency (e.g., [.01 .25]), then TR.
%                 In all cases, after the frequency value, you need to provide TR in sec.
%        - *Example:* 'hpf', .01, TR
%                     'bpf', [.01 .25], TR
%
%   **regressors:**
%        - This option will run additional GLM with the regressors specified
%          by the experimenter. The number of rows should be the number of
%          images, and the number of columns are the number of regressors. 
%        - with extract_roi, the function will extract roi values from beta
%          images as well.
%        - Example: 'regressors', X (you can build your own design matrix,
%                                   e.g., hrf convolution using onset and 
%                                   duration, using onsets2fmridesign.m)
%
%   **extract_roi:**
%        This option will extract data from ROIs specified. This
%        option should be followed by one or more masks.
%        For one mask (potentially multiple ROIs, enter a char array with the mask name.
%        For multiple masks (1 or more), enter in a cell array of mask names.
%        You can specify methods with 'roi_methods' option.
%        - 'average_over' (default): calculate averaged value across the ROIs.
%        - 'pattern_expression': calculate dot-products between
%                                pattern mask and data
%        - 'unique_mask_values' (default): will divide a mask into
%                                          multiple regions that have different discrete values.
%        - 'contiguous_regions': will divide a mask into multiple
%                                contiguous regions.
%        - 'whole': will do average_over or pattern_expression across
%                   all the voxels within the mask.
%        - *Example:* 'extract_roi', mask, 'contiguous_regions'
%                     'extract_roi', mask, 'pattern_expression'
%
%   **no_preproc:**
%        If you want to skip the preprocessing part, and want to
%        extract ROI values only, you can use this option.
%
%   **no_plots:**
%        To suppress plots, enter 'no_plots'
%
% :Outputs:
%
%  **preprocessed_dat:**
%        fmri_data object after removing nuisance variables and
%        filtering temporal confounds.
%
%  **roi_val:**
%        returns values extracted from ROIs in cell arrays (if there are many different ROIs).
%        Each cell will have roi_val.dat, roi_val.mask_name, and roi_val.methods.
%
%  **maskdat:**
%        returns resampled roi (or pattern) masks 
%
%  **beta_dat:**
%        returns beta images from the additional GLM (see 'regressors' option)
%
%
% :Examples:
% ::
%
%     roi_masks = which('weights_NSF_grouppred_cvpcr.img');
%     [preprocessed_dat, roi_val] = canlab_connectivity_preproc(dat, 'vw', 
%               'bpf', [.008 .25], TR, 'extract_roi', roi_masks, 'pattern_expression');


% ..
%    PROGRAMMER'S NOTE
%    05/19/15 fixed a bug related to conn_filter
%    08/24/18 implemented simultaneous pass filter + regression
%    01/05/19: Tor and Marianne: Pass filter all covariates (even
%               spikes., Consistent with principles in Lindquist 2019
%               "sequential filtering" paper. See more notes below.
%    08/01/19: Wani reversed Tor and Marianne's revision above 
%               after dicussion with Tor and others. Spikes should be removed
%               before "low-pass filtering" the nuisance covariates, and
%               added to the covariates later (see Ciric et al., 2018)
%               Lindquist et al., 2019 did not address the issue related to 
%               the use of spike regressors combined with low pass
%               filtering. This issue should be further examined later. 
%    08/01/19: Wani added an option for additional GLM
% ..

additional_R = [];
remove_vent_wm = false;
use_canonical_masks = true;
do_filter = false;
do_extract_roi = false;
do_preproc = true;
do_windsorize = false;
do_linear = false;
do_plots = true;
regressors = [];
do_additional_regress = false;

z = '===============================================================';
zz = '---------------------------------------------------------------';
fprintf('%s\n', z);

% parsing varargin
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'additional_nuisance', 'additional_R'}
                additional_R = varargin{i+1};
                
            case {'vw', 'ventricle_whitematter', 'remove_vw'}
                remove_vent_wm = true;
                
                % filter
            case {'hpf', 'lpf', 'bpf'}
                do_filter = true; 
                TR = varargin{i+2};
                
                if any(strcmp(varargin, 'hpf')), filter(1) = varargin{i+1}; filter(2)= Inf; end
                if any(strcmp(varargin, 'lpf')), filter(1) = 0; filter(2) = varargin{i+1}; end
                if any(strcmp(varargin, 'bpf')), filter = varargin{i+1}; end
                
            case {'windsorize'}
                do_windsorize = true;
                k_std = varargin{i+1};
                
                % detrend
            case {'linear_trend'}
                do_linear = true;
                
                % extract roi values
            case {'extract_roi', 'roi', 'mask'}
                do_extract_roi = true;
                mask = varargin{i+1};
                
            case {'no_preproc'}
                do_preproc = false;
                
            case {'no_plots'}
                do_plots = false;
                
            case {'regressors'}
                regressors = varargin{i+1};
                do_additional_regress = true;
                regressors(:,all(regressors==1)) = []; % remove an intercept if any
                nreg = size(regressors,2);
        end
    end
end

%% Remove nuisance

if do_preproc
    
    tic;
    
    nobs = size(dat.dat, 2);
    if do_additional_regress
        if size(regressors,1) ~= nobs
            error('The number of rows of regressors is different from the number of images from dat. Please check.');
        end
    end
    
    fprintf('Canlab_connectivity_preproc: starting preprocessing ... \n');
    comb_R = [regressors dat.covariates additional_R]; % if all these are empty, comb_R should still be empty
    
    if do_additional_regress
        beta_dat = dat;
        beta_dat.dat = [];
    end
    
    dat.dat = double(dat.dat);  % enforce double to avoid data format problems
    dat.history{end+1} = 'Processed with canlab_connectivity_preproc';
    
    %% white matter, ventricle
    
    if do_plots
        create_figure('plot_qc',2,4);
        set(gcf, 'position', get(0,'screensize')); drawnow;
        plot_qc(dat,1,0,[]); drawnow;
    end
    
    if remove_vent_wm
        
        fprintf('::: removing ventricle & white matter signal ... \n');
        use_pca_comp_vw = true; % default
        
        if any(strcmp(varargin, 'raw')), use_pca_comp_vw = false; end
        
        if any(strcmp(varargin, 'datdir'))
            use_canonical_masks = false;
            subject_dir = varargin{find(strcmp(varargin, 'datdir'))+1}; 
        end
        
        if use_canonical_masks 
            
            [value, components] = extract_gray_white_csf(dat, 'masks', ...
                {'gray_matter_mask.nii', 'canonical_white_matter_thrp5_ero1.nii', ...
                'canonical_ventricles_thrp5_ero1.nii'});
            wm_nuisance = value(:,2:3); % white, CSF
            wm_nuisance_comps = [components{2} components{3}]; % 5 components for white, CSF
        
        else
            
            if ~exist(fullfile(subject_dir, 'Structural/SPGR/white_matter.img'), 'file')...
                    && ~exist(fullfile(subject_dir, 'Structural/SPGR/ventricles.img'), 'file')
                wm_mask = filenames(fullfile(subject_dir, 'Structural/SPGR/wc2*.nii'), 'char', 'absolute');
                gm_mask = filenames(fullfile(subject_dir, 'Structural/SPGR/wc1*.nii'), 'char', 'absolute');
                canlab_create_wm_ventricle_masks(wm_mask, gm_mask, 'wm_thr', .9, 'vent_thr', .9);
            end
            
            [wm_nuisance, wm_nuisance_comps] = canlab_extract_ventricle_wm_timeseries(fullfile(subject_dir, 'Structural/SPGR'), dat, 'noplot');
        
        end
        
        if use_pca_comp_vw
            
            comb_R = [comb_R scale(wm_nuisance_comps)];
            
        elseif ~use_pca_comp_vw
            
            comb_R = [comb_R scale(wm_nuisance)];
            
        end
        
        dat.history{end+1} = 'Extracting white matter and ventricle signals';
        
    end
    
    %% linear trend
    if do_linear
        
        comb_R(:,end+1) = scale((1:nobs)');
        dat.history{end+1} = '(added linear detrending to nuisance covs)';
        
    end
    
    % (optional) global signal % not implemented
    
    if do_filter

            fprintf('::: temporal filtering brain data... \n');
            % conn_filter does mean-center by default.
            % We need mean back in
            
            % Pass-filter brain data
            % -------------------------------------------
            dat.dat = conn_filter(TR, filter, dat.dat', 'full')' ...
                + repmat(mean(dat.dat,2), 1, size(dat.dat,2)); 
            
            dat.history{end+1} = 'Done temporal filtering brain data';
            
            if do_plots
                plot_qc(dat,2,1,[]); drawnow;
            end
    end
    
    %% get residuals
    
    if ~isempty(comb_R)
        % Do regression/residualization
        
        if size(comb_R, 1) ~= nobs
            error('Nuisance covariates and data matrix have different nunmbers of cases.');
        end
        
        %% temporal filter (low and high): low pass, high pass, bandpass
        % Simultaneous temporal filtering and regression
        
        if do_filter
            
            fprintf('::: temporal filtering covariates ... \n');
            % conn_filter does mean-center by default.
            % We need mean back in
            
            comb_R = filter_cov_wo_spikes(comb_R, TR, filter);
            
            dat.history{end+1} = 'Done temporal filtering covariates';
            
        end
        
        % add intercept
        X = [comb_R ones(nobs,1)];
        
        if do_additional_regress
            
            % regression
            beta_dat.dat = data.dat * pinv(X)';
            beta_dat.dat = beta_dat.dat(:,1:nreg);
            beta_dat.covariates = X;
            X(:,1:nreg) = [];
        end
        
        % residualize
        dat.dat = dat.dat - dat.dat * pinv(X)' * X'; % this is more efficient than (dat.dat' - X*pinv(X)*dat.dat')'
        dat.history{end+1} = 'Regressed out nuisance covariates';
        
        if do_plots
            clim = plot_qc(dat,3,1,[]); drawnow;
        end
        
        
    else  % no regression
        
    end
    
    %% Windsorize entire data matrix to k STD
    if do_windsorize
        dat = windsorize(dat, k_std); % history added within the function
        if do_plots
            plot_qc(dat,4,0,clim); drawnow;
        end
    end
    
    t = toc;
    fprintf('Done: data preprocessing in %4.1d sec.\n\n', t);
    fprintf('%s\n\n', zz);
end

%% get ROI values

if do_extract_roi
    tic;
    fprintf('Canlab_preproc: extracting ROI signal ... \n');
    
    % default
    do_pattern = false;
    do_whole = false;
    unique_mask_values = true;
    contiguous_regions = false;
    
    if any(strcmp(varargin, 'pattern_expression')), do_pattern = true; unique_mask_values = false; end
    if any(strcmp(varargin, 'whole')), do_whole = true; unique_mask_values = false; end
    if any(strcmp(varargin, 'contiguous_regions')), contiguous_regions = true; unique_mask_values = false; end
    
    if ~iscell(mask)
        mask_temp = mask; clear mask;
        mask = cell(size(mask_temp,1),1); % preallocation
        for i = 1:size(mask_temp,1)
            mask{i} = mask_temp(i,:);
        end
    end
    
    [roi_val, maskdat] = extract_roi_sub(dat, mask, do_pattern, do_whole, unique_mask_values, contiguous_regions);
    
    if do_additional_regress
        beta_roi_val = extract_roi_sub(beta_dat, mask, do_pattern, do_whole, unique_mask_values, contiguous_regions);
    end
    
    
    % output: roi_val, maskdat
    % input: mask, dat, do_pattern, do_whole, unique_mask_values, contiguous_regions
    
    t2 = toc;
    fprintf('Done in %4.1d sec.\n', t2);
    fprintf('%s\n', zz);
    
    if do_preproc
        fprintf('Total time: 4.1%d sec.\n', t+t2);
        fprintf('%s\n', z);
    end
else
    roi_val = {};
end

end


%%  ---------------------------------------------------
%    Sub-functions
%   ---------------------------------------------------

function [y,fy]=conn_filter(rt, filter, x, option)

% from conn toolbox
% http://www.nitrc.org/projects/conn/

if nargin < 4, option='full'; end

fy=fft(x,[],1);
f=(0:size(x,1)-1);
f=min(f,size(x,1)-f);

switch(lower(option))
    case 'full'
        
        idx=find(f<filter(1)*(rt*size(x,1))|f>=filter(2)*(rt*size(x,1)));
        %idx=idx(idx>1);
        fy(idx,:)=0;
        y=real(ifft(fy,[],1))*2*size(x,1)*(min(.5,filter(2)*rt)-max(0,filter(1)*rt))/max(1,size(x,1)-numel(idx));
        
    case 'partial'
        
        idx=find(f>=filter(1)*(rt*size(x,1))&f<filter(2)*(rt*size(x,1)));
        %if ~any(idx==1), idx=[1,idx]; end
        y=fy(idx,:);
        
end

end

function obj = windsorize(obj, k)

% Windsorize entire data matrix to k x STD

m = sum(sum(obj.dat)) ./ numel(obj.dat); % save memory space

s = k .* std(obj.dat(:));

whbad = obj.dat < m - s;
obj.dat(whbad) = m - s;
whbad2 = obj.dat > m + s;
obj.dat(whbad2) = m + s;

nbad = sum(whbad(:)) + sum(whbad2(:));
percbad = 100 .* nbad ./ numel(obj.dat);

obj.history{end+1} = sprintf('Windsorized data matrix to %d STD; adjusted %3.0f values, %3.1f%% of values', k, nbad, percbad);
disp(obj.history{end});
end

function clim = plot_qc(dat, n, use7sd, clim)

% dat
% n: steps
% use7sd: use -6sd~+6sd as a range
% clim: specify color limit

titles = {'RAW:before Conn preproc', 'Pass filter', 'Residualize (Nuisance, VentWM, linear trend, if selected)', 'Windsorize'};

subplot(2, 4, n);
plot(mean(dat.dat)); 
set(gca, 'xlim', [0 size(dat.dat,2)]);
xlabel('timeseries');
ylabel('Global mean');
title(titles{n});

subplot(2, 4, n+4);

if use7sd
    m = nanmean(dat.dat(:));
    sd = nanstd(dat.dat(:));
    clim = [m-7*sd m+7*sd];    
end

if isempty(clim)
    imagesc(dat.dat);
else
    imagesc(dat.dat, clim);
end

colormap(gray);
set(gca, 'xlim', [0 size(dat.dat,2)]);
xlabel('timeseries');
ylabel('voxels');
colorbar('southoutside');

end

function R = filter_cov_wo_spikes(R, TR, filter)

% Filter the covariates
% Eliminate bias from sequential filtering operations
% -------------------------------------------
% 1/15/2019: Tor and Marianne: Pass filter all covariates (even
% spikes), because filter transforms data in a way that matches
% what is done to spike regressors. So even if spike regressors
% "look funny", they adjust for the effects of "spikes" on the
% time series data appropriately, accounting for filtering.
%
% 8/1/2019: Wani reversed this to the previous version after
% dicussion with Tor and others. Spikes should be removed
% before "low-pass filtering" the nuisance covariates, and
% added to the covariates later (see Ciric et al., 2018)

% pass-filter covs (without outlier indices)
outlier_idx = sum(R==1)==1;
outliers = R(:,outlier_idx);
R(:,outlier_idx) = [];

R = conn_filter(TR, filter, R, 'full') ...
    + repmat(mean(R), size(R,1), 1);

% add outlier indices back in
R = [R outliers];

end

function [roi_val, maskdat] = extract_roi_sub(dat, mask, do_pattern, do_whole, unique_mask_values, contiguous_regions)

% subfunction for extract rois

% copy over apply_mask
maskdat = cell(numel(mask),1);
for i = 1:numel(mask)
    maskdat{i} = fmri_data(mask{i});
    isdiff = compare_space(dat, maskdat{i});
    
    if isdiff == 1 || isdiff == 2 % diff space, not just diff voxels
        maskdat{i} = resample_space(maskdat{i}, dat, 'nearest');
        
        if length(maskdat{i}.removed_voxels) == maskdat{i}.volInfo.nvox
            disp('Warning: resample_space returned illegal length for removed voxels. Fixing...');
            maskdat{i}.removed_voxels = maskdat{i}.removed_voxels(maskdat{i}.volInfo.wh_inmask);
        end
    end
end

roi_val = cell(numel(mask),1); % preallocation

for i = 1:numel(mask)
    if do_whole
        % single, whole-brain pattern
        if ~do_pattern
            
            roi_obj = apply_mask(dat, maskdat{i});
            roi_obj = remove_empty(roi_obj);
            roi_val{i}.dat = nanmean(roi_obj.dat)';
            roi_val{i}.mask_name = mask{i};
            roi_val{i}.methods = 'averaged_over using a whole mask';
            
        elseif do_pattern
            
            roi_val{i}.dat = apply_mask(dat, maskdat{i}, 'pattern_expression', 'ignore_missing');
            roi_val{i}.mask_name = mask{i};
            roi_val{i}.methods = 'pattern expression using a whole mask';
            
        end
    else
        % separate patterns for local regions
        
        if ~do_pattern
            if unique_mask_values
                
                roi_obj = extract_roi_averages(dat, maskdat{i});
                roi_val{i}.dat = cat(2, roi_obj(:).dat);
                
                % extract_roi_averages does not preserve missing cases: remove
                if any(dat.removed_images)
                    roi_val{i}.dat = roi_val{i}.dat(~dat.removed_images, :);
                end
                
                roi_val{i}.mask_name = mask{i};
                roi_val{i}.methods = 'averaged_over using unique_mask_values';
                
            elseif contiguous_regions
                
                roi_obj = extract_roi_averages(dat, maskdat{i}, 'contiguous_regions');
                roi_val{i}.dat = cat(2, roi_obj(:).dat);
                
                % extract_roi_averages does not preserve missing cases: remove
                if any(dat.removed_images)
                    roi_val{i}.dat = roi_val{i}.dat(~dat.removed_images, :);
                end
                
                roi_val{i}.mask_name = mask{i};
                roi_val{i}.methods = 'averaged_over using contiguous_regions';
                
            else
                warning('Unknown ROI extraction option. Not extracting ROI data...');
            end
        else
            % local patterns within each region
            
            roi_obj = extract_roi_averages(dat, maskdat{i}, 'pattern_expression', 'contiguous_regions');
            roi_val{i}.dat = cat(2, roi_obj(:).dat);
            
            % extract_roi_averages does not preserve missing cases: remove
            if any(dat.removed_images)
                roi_val{i}.dat = roi_val{i}.dat(~dat.removed_images, :);
            end
            
            roi_val{i}.mask_name = mask{i};
            roi_val{i}.methods = 'pattern expression using contiguous_regions';
        end
    end
end

end
