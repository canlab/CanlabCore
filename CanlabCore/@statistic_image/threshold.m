function stats_image_obj = threshold(stats_image_obj, input_threshold, thresh_type, varargin)
% Threshold statistic_image object based on statistical threshold values.
%
% :Usage:
% ::
%
%    stats_image_obj = threshold(stats_image_obj, pvalthreshold or other thresh, thresh_type, ['k', extent_thresh OR other optional inputs])
%
% This is a method for an statistic_image object.
% Thresholding is reversible.
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2015 Tor Wager
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
%   **stats_image_obj:**
%        statistic_image object
%
%   **input_threshold:**
%        [pvalthreshold or other thresh]
%        A numeric value corresponding to the threshold desired.
%        Either a p-value or a range of raw values, depending on the threshold
%        type.
%
%
%   **thresh_type:**
%        Threshold type which can be one of:
%          - 'fdr' : FDR-correct based on p-values already stored in image .p field
%          - 'fwe' : FWE-correct; not implemented
%          - 'bfr' : Bonferroni correction (FWE).
%          - 'unc' : Uncorrected p-value threshold: p-value, e.g., .05 or .001
%          - 'extent', 'cluster_extent' : Cluster extent correction with GRF at p < .05 corrected, 
%                   primary threshold determined by input_threshold
%          - 'raw-between' : threshold raw image values; save those > input_threshold(1) and < input_threshold(2)
%          - 'raw-outside' : threshold raw image values; save those < input_threshold(1) or > input_threshold(2)
%
% :Optional Inputs:
%
%   **k:**
%        Followed by cluster extent in voxels: extent-based thresholding of any of the above
%
%   **noverbose:**
%        Suppress verbose output
%
%   **mask:**
%        Followed by name of mask or fmri_mask_image object
%          - this will affect corrected significance levels
%
% :Outputs:
%
%  **stats_image_obj:**
%        thresholded statistic_image object
%
% :Example:
% ::
%
%    % Retain sig pos or neg results at p < .001 uncorrected, cluster extent >= 100 voxels
%    obj = threshold(obj, .001, 'unc', 'k', 100)
%
%    % Retain sig pos or neg results at q < .05 FDR, cluster extent >= 10 voxels
%    obj = threshold(obj, .05, 'fdr', 'k', 10)
%
%    % Retain voxels with absolute statistic/data value > 3
%    obj = threshold(obj, [-3 3], 'raw-outside')
%
%    dat = threshold(dat, 0.001, 'unc', 'k', 35, 'mask', which('scalped_avg152T1_graymatter_smoothed.img'));
%    dat = threshold(dat, 0.001, 'unc', 'k', 35, 'mask', maskobj);
%
% Here is a complete example using simulated data.  First, we load a
% canonical network and generate true-signal areas. Second, do a t-test on the simulated images.
% Then, load two different masks and perform cluster extent correction within one of
% the masks.
% imgs = load_image_set('bucknerlab');
% true_mask = get_wh_image(imgs, 2); % select a network for true signal
% patient_obj = sim_data(fmri_data, 'n', N, 'd', 1, 'smoothness', 10, 'true_region_mask', true_mask);
% mask1 = fmri_data(which('gray_matter_mask.img'), 'noverbose');
% mask2 = fmri_data(which('brainmask.nii'), 'noverbose');
% patient_t = ttest(patient_obj);
% patient_t = threshold(patient_t, .01, 'cluster_extent', patient_obj, 'mask', mask2);
% 
% 
% :See also:
% image_vector.threshold, statistic_image.multi_threshold
%
% ..
%    Programmers' notes:
%    Created by Tor Wager, Dec 2010
%    Tor: Updated documentation, July 2015
%
%    7/19/2013: added 'mask' and 'bfr'(bonferroni) option by Wani Woo.
%    With the 'mask' option, you can define a space for the multiple comparison
%    correction.
% ..

k = 1;
doverbose = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'k'
                k = varargin{i + 1};
                
                if isempty(stats_image_obj.volInfo)
                    error('You must add a volInfo structure to the statistic image object to do extent-based thresholding');
                end
                
            case 'mask' % 7/19/2013: added by Wani Woo.
                
                maskinput = varargin{i + 1}; varargin{i + 1} = 0;
                switch class(maskinput)
                    case 'char' % string file name
                        % maskobj = fmri_mask_image(maskinput);
                        maskobj = fmri_data(maskinput, 'noverbose');
                        
                    case {'fmri_mask_image', 'fmri_data', 'image_vector'}
                        maskobj = maskinput;
                        
                    otherwise
                        error('region class constructor: unknown mask input type.')
                end
                clear maskinput
                
                % Tor added: check and/or perform resampling if needed
                maskobj = resample_space(maskobj, stats_image_obj);
                
                maskobj = replace_empty(maskobj);
                stats_image_obj = replace_empty(stats_image_obj);
                
                % Check
                if length(maskobj.dat) ~= length(stats_image_obj.dat)
                    if sum(sum(maskobj.volInfo.mat ~= stats_image_obj.volInfo.mat))
                        mask_mapdat = scn_map_image(maskobj, stats_image_obj);
                        maskdat = logical(mask_mapdat(stats_image_obj.volInfo.wh_inmask));
                    else
                        error('VolInfo of mask and data are same, but data lengths are different.');
                    end
                else
                    if sum(sum(maskobj.volInfo.mat ~= stats_image_obj.volInfo.mat))
                        error('VolInfo of mask and data are diffferent, but data lengths are same.');
                    else
                        maskdat = logical(maskobj.dat);
                    end
                end
                
                stats_image_obj = apply_mask(stats_image_obj, maskobj);
                
%                 stats_image_obj = replace_empty(stats_image_obj);
%                 stats_image_obj.p = stats_image_obj.p(maskdat);
%                 stats_image_obj.dat = stats_image_obj.dat(maskdat);
%                 if ~isempty(stats_image_obj.sig)
%                     stats_image_obj.sig = stats_image_obj.sig(maskdat);
%                 end
%                 stats_image_obj.removed_voxels = ~maskdat;
                
            case 'noverbose'
                doverbose = 0;
              
                %             case {'extent', 'cluster_extent'}
                %
                %                 resid_image_obj = varargin{i + 1}; varargin{i + 1} = [];
                
            otherwise
                error('Illegal argument keyword.');
        end
    end
end

% Extent thresholding
switch thresh_type
    case {'extent', 'cluster_extent'}
        
        resid_image_obj = varargin{1};
        
    otherwise
end


n = size(stats_image_obj.dat, 2); % images are columns here

for i = 1:n
    
    % Determine threshold
    % ---------------------------------------------------
    
    switch lower(thresh_type)
        case 'fdr'
            % fix zero P-vals, which can occur with some nonparam methods
            % or large t-values
            wh = stats_image_obj.p == 0 & abs(stats_image_obj.dat) > 1000*eps;
            stats_image_obj.p(wh) = 1000*eps;
            
            thrval = FDR(stats_image_obj.p(:, i), input_threshold);
            if isempty(thrval)
                thresh(i) = 0;
            else
                thresh(i) = thrval;
            end
            fprintf('Image %3.0f FDR q < %3.3f threshold is %3.6f\n', i, input_threshold, thresh(i));
            stats_image_obj.threshold(i) = thresh(i);
            
        case 'fwe'
            error('not implemented yet.')
            
        case {'bfr', 'bonferroni'} % 7/19/2013: added by Wani Woo.
            thresh(i) = input_threshold/size(stats_image_obj.p,1);
            stats_image_obj.threshold(i) = thresh(i);
            
        case {'unc', 'uncorrected'}
            thresh(i) = input_threshold;
            stats_image_obj.threshold(i) = thresh(i);
            
        case {'raw-between', 'raw-outside'}
            thresh = input_threshold;
            
        case {'extent', 'cluster_extent'}
            
            % Cluster extent thresholding with SPM's GRF
            % At primary threshold  = input_threshold, p < .05 corrected
            % Recommended mask: which('gray_matter_mask.img')
            
            % cl_ext_spm_grf takes image file names as inputs, so create
            % files.
            % Temporary mask - if mask was entered, we have already masked
            % this object above, so use all voxels here as mask (NOT just
            % ones above previous threshold)
           
            extent_mask = stats_image_obj;
%             n_in_mask = size(extent_mask.dat, 1);
%             extent_mask.dat = ones(n_in_mask, 1);
%             extent_mask.sig = true(n_in_mask, 1);
%             extent_mask.ste = ones(n_in_mask, 1);
%             extent_mask.N = ones(n_in_mask, 1);
%             extent_mask.removed_voxels = [];
            extent_mask.fullpath = fullfile(pwd, 'tmp_mask_img_remove_me.nii');
            write(extent_mask);

            % Write images to disk - we will remove after thresholddetermination
            resid_image_obj.image_names = 'tmp_resid_imgs_remove_me.nii';
            resid_image_obj.fullpath = fullfile(pwd, resid_image_obj.image_names);
            write(resid_image_obj);

            % Get extent threshold
            [k, fwhm] = cl_ext_spm_grf(.05, input_threshold, resid_image_obj.fullpath, extent_mask.fullpath);
            
            thresh(i) = input_threshold; % for later
            
            if doverbose
                fprintf('Cluster extent correction: Estimated smoothness FWHM is %3.2f x %3.2f x %3.2f\nExtent threshold k = %3.0f\n', fwhm, k);
            end
            
            if n > 1
                disp('***Warning!!! Cluster extent correction requires a set of residual images for each specific effect.')
                disp('You are thresholding multiple images with one set of residual images, which may be invalid.')
                disp('Separating images in the object you are thresholding and entering separate resid_image_objects is a solution.');
            end
            
            delete(resid_image_obj.fullpath);
            delete(extent_mask.fullpath);

        otherwise
            error('Unknown threshold type. Enter fdr, fwe, or uncorrected, or raw-between or raw-outside')
    end
    
    %     % Re-set .dat depending on image type (needs to be updated for other
    %     % image types)
    %     % ---------------------------------------------------
    %     switch
    
    switch lower(thresh_type)
        
        case {'fdr', 'fwe', 'unc', 'uncorrected', 'bfr', 'bonferroni', 'extent', 'cluster_extent'}
            stats_image_obj.sig(:, i) = stats_image_obj.p(:, i) < thresh(i);
            
            % Apply size threshold
            % --------------------------------------
            stats_image_obj = replace_empty(stats_image_obj);
            
            if k > 1
                
                % reparse contig vox to consider only those above threshold

                % stats_image_obj = reparse_contiguous(stats_image_obj);
%                 stats_image_obj = replace_empty(stats_image_obj);
                % iimg_cluster_extent not always working right - 7/27/2018
                % - under some strange cases i can't quite understand (tor)
                
                %stats_image_obj.sig(:, i) = logical(iimg_cluster_extent(double(stats_image_obj.sig(:, i)), stats_image_obj.volInfo, k));
            
                % apply the .sig field, re-derive cluster indices in sig
                % voxels, and include only clusters with enough voxels.
                stats_image_obj = alternate_cluster_extent_in_sig_field(stats_image_obj, i, k);
                
            end
            
            if doverbose
                
                [n_vox_in_cluster, indic, cluster_indx_vals] = get_cluster_sizes(stats_image_obj);
                wh = n_vox_in_cluster > 0;
                
                fprintf('\nImage %3.0f\n%3.0f contig. clusters, sizes %3.0f to %3.0f\n', ...
                    i, sum(wh), ...
                    min(n_vox_in_cluster(wh)), max(n_vox_in_cluster(wh)));
                
                fprintf('Positive effect: %3.0f voxels, min p-value: %3.8f\n', ...
                    sum(stats_image_obj.dat(:, i) > 0 & stats_image_obj.sig(:, i)), ...
                    min(stats_image_obj.p(stats_image_obj.dat(:, i) > 0, i)));
                
                fprintf('Negative effect: %3.0f voxels, min p-value: %3.8f\n', ...
                    sum(stats_image_obj.dat(:, i) < 0 & stats_image_obj.sig(:, i)), ...
                    min(stats_image_obj.p(stats_image_obj.dat(:, i) < 0, i)));
                
            end
            
        case 'raw-between'
            stats_image_obj.sig(:, i) = stats_image_obj.dat(:, i) > thresh(1) & stats_image_obj.dat(:, i) < thresh(2);
            
            % Apply size threshold
            % --------------------------------------
            if k > 1
                stats_image_obj.sig(:, i) = logical(iimg_cluster_extent(double(stats_image_obj.sig(:, i)), stats_image_obj.volInfo, k));
            end
            
            if doverbose
                fprintf('Keeping vals between %3.3f and %3.3f: %3.0f voxels in .sig\n', thresh(1), thresh(2), sum(stats_image_obj.sig(:, i)));
            end
            
        case 'raw-outside'
            stats_image_obj.sig(:, i) = stats_image_obj.dat(:, i) < thresh(1) | stats_image_obj.dat(:, i) > thresh(2);
            
            % Apply size threshold
            % --------------------------------------
            if k > 1
                stats_image_obj.sig(:, i) = logical(iimg_cluster_extent(double(stats_image_obj.sig(:, i)), stats_image_obj.volInfo, k));
            end
            
            if doverbose
                fprintf('Keeping vals outside of %3.3f to %3.3f: %3.0f voxels in .sig\n', thresh(1), thresh(2), sum(stats_image_obj.sig(:, i)));
            end
            
        otherwise
            error('Unknown threshold type. Enter fdr, fwe, or uncorrected')
    end
    
    % Replace .dat field with zero for non-sig voxels, so that we can
    % create regions, etc.  .t field and .p field, etc., still have
    % original info.
    % But doing this here breaks the ability to re-threshold, so now done
    % in region.m
    %stats_image_obj.dat(~stats_image_obj.sig(:, i), i) = 0;
    
    
end  % image loop

if isempty(stats_image_obj.volInfo)
    if doverbose
        disp('Warning! volInfo not defined for stats image object.  regions cannot be formed correctly from this image.')
    end
else
    
    
    % Update volInfo, for region-defining purposes.
    % Note that this will only work for the FIRST image in the set
    
    stats_image_obj = reparse_contiguous(stats_image_obj);
    
end

end % function




function stats_image_obj = alternate_cluster_extent_in_sig_field(stats_image_obj, i, k)


stats_image_obj = reparse_contiguous(stats_image_obj);  % clusters only for significant vox in .sig, or 0.

stats_image_obj = replace_empty(stats_image_obj);       % sig is now all in-mask. also reparse_contiguous may remove empties...

% mysig = logical(stats_image_obj.sig(:, i));             % significant vox vector to prune (remove vox where clusters too small)
new_sig = false(size(stats_image_obj.sig(:, i)));

new_cluster = zeros(size(stats_image_obj.volInfo.cluster)); 

obj_tmp = get_wh_image(stats_image_obj, i);             % Use sig vox for this region only in reparsing contig and counting...
[n_vox_in_cluster, indic, cluster_indx_vals] = get_cluster_sizes(obj_tmp);

wh_include = n_vox_in_cluster >= k;

% build up new sig, including only voxels where cluster index is not 0
% after masking with original mask in .sig 
for j = 1:length(wh_include)
    
    if wh_include(j)
        
        new_sig(indic(:, j), 1) = true;
        
        new_cluster(indic(:, j), 1) = cluster_indx_vals(j); % this does not renumber clusters. leaves them intact.
    end
    
end

stats_image_obj.sig(:, i) = new_sig;
stats_image_obj.volInfo.cluster = new_cluster;

end % function




function [n_vox_in_cluster, indic, cluster_indx_vals] = get_cluster_sizes(stats_image_obj)
%
%
% n_vox_in_cluster: Counts of voxels in each contiguous cluster indexed 1 to k, as indexed by stats_image_obj.volInfo.cluster
%                   If index values are missing, counts can be zero.
% 
% indic:            Voxels x index values logical matrix showing which voxels belong to which cluster
%
% cluster_indx_vals: Vector of integers with index values for clusters. May not be consecutive if some indices are missing (empty clusters)
 
% Count contiguous clusters
% ---------------------------------------------------------------
stats_image_obj = reparse_contiguous(stats_image_obj);            % clusters only for significant vox in .sig, or 0.
                                                                  % re-do because this may be called as a stand-alone function

cindx = stats_image_obj.volInfo.cluster;                          % contig cluster index, all in-mask, zero for not-.sig

[indic, cluster_indx_vals] = condf2indic(cindx, 'integers');      % omit 0, list of all integers for cluster IDs

indic(isnan(indic)) = 0;

indic = logical(indic);

n_vox_in_cluster = double(sum(indic));
% ---------------------------------------------------------------

end

