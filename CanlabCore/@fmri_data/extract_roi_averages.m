function [cl, clroimean, clpattern] = extract_roi_averages(obj, mask_image, varargin)
% This fmri_data method a extracts and averages data stored in an fmri_data object
% from a set of ROIs defined in a mask.
%
% If no mask_image is entered, it uses the mask defined with the fmri_data object as a default.
%
% If mask_image is a new image file name, this method:
%   1. Defines an fmri_mask_image object using create_fmri_mask
%   2. Maps to the space in fmri_data object using resample_to_image_space
%
% Regions to average over can be either regions of contiguous voxels
% bounded by voxels with values of 0 or NaN, which are considered non-data
% values, or regions defined by unique integer codes in the mask image
% (i.e., for atlas images with unique codes for each defined region.)
%
% Mask/Atlas image does NOT have to be in the same space as the images to
% extract from.  It will be remapped/resliced.
% NOTE: Mask is *reloaded* from original data if space is remapped, and you
% cannot use manual thresholding of the mask. This is a feature of the
% map_to_image_space method and scn_map_image
%
% Extracted data is returned in single data format.
%
% :Usage:
% ::
%
%     [[cl, clroimean, clpattern] = extract_roi_averages(fmri_data obj, [mask_image], [average_over])
%
% :Inputs:
%   1. char array of strings containing 4D image file names (data extracted from these)
%   2. mask_image to extract from.
%
% :Optional inputs:
%   **how to average**
%       Default = 'unique_mask_values' to average over unique integer codes in the mask image
%       bounded by voxels of 0 or NaN (non-data values)
%       (i.e., for atlas images with unique codes for each defined region)
%       Alt. option = 'contiguous_regions' to average over contiguous voxels
%
%   **pattern_expression**
%       Use values in mask images to get weighted average within each
%       region, rather than simple average.  See also apply_mask with
%       'pattern_expression' option.
%
%       Optional outputs (varargout):
%       [cl, cl_roimean, cl_roipattern] = ...
%       roimean: pattern expression is average over ROI (unit vector)
%       roipattern: pattern expression is dot product of activity and mean-centered pattern weights
%
%   **nonorm**
%       Turn off L1 norm in pattern expression.
%
%   **cosine_similarity**
%       Use cosine similarity metric for pattern expression instead of dot product
%       Passed in to canlab_pattern_similarity
%
%   **correlation**
%       Use correlation metric for pattern expression instead of dot product
%       Passed in to canlab_pattern_similarity
%
% :Examples:
% ::
%
%     imgs_to_extract_from = filenames('w*.nii','char');
%     mask_image = which('anat_lbpa_thal.img');
%     [cl, clroimean, clpattern] = extract_image_data(imgs_to_extract_from, mask_image);
%
%     region_obj = extract_roi_averages(data_obj, mask_char_name, 'pattern_expression', 'contiguous_regions');
%
% :Related functions:
%    For an non-object-oriented alternative, see extract_image_data.m
%    apply_parcellation.m
%
% ..
%     Notes:
%     cl(i).dat gives you the pattern expression values for cluster i.
%
%     This function LOSES removed image data - you must re-remove if you have
%     removed images!
%
%     Modified June 11, 2013 by Tor
%       - use resample_space instead of resample_to_image_space
%
%     Modified Dec 1, 2014 by Wani
%       - moved up the part of parsing optional inputs because resample_space
%         causes a problem for the unique_mask_values option
%       - For resample_space, the 'nearest' option should be used when the
%         "unique_mask_values" option is used.
%
%     Modified Oct 2015 by Tor]
%       - Clarified options, empty cl error check, changed varargout behavior
% ..
% Programmers' notes:
% This function is different from fmri_data.extract_roi_averages
% Better to have only one function of record in the future...
% Note: 
% cl = extract_roi_averages(imgs, atlas_obj{1});
% accomplishes the same task as apply_parcellation, returns slightly different values due to interpolation.  
% Also to-do: Uniform code with match_spaces in apply_parcellation.

pattern_norm = 1; % for pattern expression -- default is norm pattern weights
doverbose = 1;

% Now handled by canlab_pattern_similarity
%similarity_metric = 'dotproduct';

% ---------------------------------
% define region object based on choices
% also define optional inputs
% ---------------------------------

average_over = 'unique_mask_values'; %'contiguous_regions'  or 'unique_mask_values';

for varg = 1:length(varargin)
    if ischar(varargin{varg})
        switch varargin{varg}
            
            % reserved keywords
            case 'contiguous_regions', average_over = 'contiguous_regions';
            case 'unique_mask_values', average_over = 'unique_mask_values';
                
            case {'pattern_expression', 'donorm'}
                % do nothing -- ignore and use later
                
            case {'nonorm'}
                pattern_norm = 0;
                
                case {'cosine_similarity', 'correlation'}
                    % Do nothing
                % Now handled by canlab_pattern_similarity

            case {'notables'}
                % Do nothing
                % Inputs passed from other functions, not needed
                
            case 'noverbose'
                doverbose = 0;
              
            case 'verbose' % do nothing

            otherwise
                disp('fmri_data.extract_roi_averages: Illegal string value for average_over.');
                fprintf('You entered ''%s''\n Valid values are %s or %s\n', varargin{varg}, '''contiguous_regions''', '''unique_mask_values''');
                error('Exiting');
        end
    end
end


%space_defining_image = deblank(obj.fullpath(1, :));
space_defining_image = obj.mask;

% ---------------------------------
% define mask object and resample to image space
% ---------------------------------

if doverbose
    fprintf('fmri_data.extract_roi_averages: ');
end

if nargin < 2 || isempty(mask_image)
    
    mask = obj.mask;
    
    obj.dat = obj.dat';
    
else
    
    if doverbose
        fprintf('Defining mask object. ');
    end
    
    switch class(mask_image)
        
        case 'char'
            % Create a mask object
            mask = fmri_mask_image(mask_image);
            
        case 'fmri_mask_image'
            mask = mask_image;
            
        case {'image_vector', 'fmri_data', 'atlas'}
            mask = replace_empty(mask_image);
            
        case 'statistic_image'
            mask = replace_empty(mask_image);
            mask.dat = mask.dat(:, 1) .* double(mask.sig(:, 1) > 0);
            mask.sig = mask.sig(:, 1);
            mask.removed_images = mask.removed_images(1);
            %mask = reparse_contiguous(mask);
            [~, vectorized_voldata] = reconstruct_image(mask);
            mask = rebuild_volinfo_from_dat(mask, vectorized_voldata);
            mask = reparse_contiguous(mask);
            
        otherwise
            error('fmri_data.extract_roi_averages: unknown mask input type.')
            
    end
    
    if doverbose && size(mask.dat, 2) > 1, fprintf('Mask contains multiple images: Using first one.\n'); end
    
    isdiff = compare_space(mask, space_defining_image);
    
    if isdiff == 1 || isdiff == 2
        % resample it to the image space
        % * this step cannot be done recursively...could be worked on.
        if doverbose, fprintf('Resampling mask. '); end
        
        if isa(mask, 'atlas')
            mask = resample_space(mask, space_defining_image); % resample_space can handle atlas objects appropriately
            
        elseif any(strcmp(average_over, 'unique_mask_values'))
            
            % added by Wani 12/01/2014: "linear" interpolation causes a
            % problem for the unique_mask_values option.
            mask = resample_space(mask, space_defining_image, 'nearest');
        else
            mask = resample_space(mask, space_defining_image);
        end
        
        %mask = resample_to_image_space(mask, space_defining_image); % do
        %not use - requires re-loading from disk and removes any manual
        %thesholding
        
    end
    
    if any(obj.volInfo.image_indx ~= mask.volInfo.image_indx)
        % space is same, but voxel indices are not the same...
        % different masks for each, most likely...
        
        % There is a simpler solution than resampling...could work on this.
        mask = resample_space(mask, obj);
        
        
    end
    
    % this mask will have a different set of in-mask voxels
    % get only the subset in BOTH the original fmri dataset mask and
    % the new one we're applying
    
    % --------------------------------------------
    % Redefine elements of new mask to apply
    % to keep only coords in original dataset
    % --------------------------------------------
    
    is_inmask = obj.mask.volInfo.image_indx & mask.volInfo.image_indx;
    mask.volInfo.image_indx = is_inmask;
    
    orig_wh_inmask = mask.volInfo.wh_inmask; % save for later -
    
    mask.volInfo.wh_inmask = find(mask.volInfo.image_indx);
    mask.volInfo.n_inmask = length(mask.volInfo.wh_inmask);
    
    % need to limit coords in new mask to ONLY those in the orig mask as
    % well.
    wh_to_keep = is_inmask(orig_wh_inmask);
    mask.volInfo.xyzlist = mask.volInfo.xyzlist(wh_to_keep, :);
    mask.volInfo.cluster = mask.volInfo.cluster(wh_to_keep, :);
    
    
    % --------------------------------------------
    % Redefine data (temporarily; not passed out)
    % to keep only coords in new mask
    % --------------------------------------------
    
    % need to limit data in obj to ONLY voxels that are also in the new
    % mask. Index in space of original fmri_data object for which to keep:
    wh_to_keep = is_inmask(obj.mask.volInfo.wh_inmask);
    
    % eliminate out-of-mask voxels before indexing into them with new mask
    obj = replace_empty(obj);
    obj.dat = obj.dat(wh_to_keep, :)';
    
    
end

if doverbose
    fprintf('\n');
end

cl = region(mask, average_over);
cl(1).source_images = obj.fullpath;

if length(cl) == 1 && isempty(cl(1).XYZ)
    % No regions!
    warning('No regions found in mask.');
    warning('The most common reason is that the average_over flag is set incorrectly.');
    fprintf('Your choice is: %s\n')
    return
end

if nargout > 1
    [clroimean, clpattern] = deal(cl);
end

% ---------------------------------
% Now get averages by cluster
% ---------------------------------

mask = replace_empty(mask); % next line will only work if replaced; tor, 5/27/15
maskData = mask.dat(logical(mask.volInfo.wh_inmask), :);

if any(strcmp(varargin, 'pattern_expression'))  % for pexp only
    maskvals = maskData;
    if doverbose, fprintf('Applying mask weights to get pattern expression.\n'); end
    
    if pattern_norm
        if doverbose, fprintf('Normalizing within contiguous region by L1 norm.\n'); end
    else
        if doverbose, fprintf('No normalization of region weights.\n'); end
    end
else
    if doverbose, fprintf('Averaging data. '); end
end

switch average_over
    
    % Define integer codes for sets of voxels to average over.
    
    case 'unique_mask_values'
        maskData = round(maskData);
        u = unique(maskData)'; u(u == 0) = [];
        nregions = length(u);
        if doverbose
            fprintf('Averaging over unique mask values, assuming integer-valued mask: %3.0f regions\n', nregions);
        end
        
    case 'contiguous_regions'
        u = unique(mask.volInfo.cluster); u(u == 0) = [];
        maskData = mask.volInfo.cluster;
        
        
    case 'none'
        cl = [];
        return
        
        
    otherwise
        error('Illegal value for average_over.  See help for this function.');
end


% Now get the average activity in each region for the defined regions

nregions = length(u);

if nregions ~= length(cl)
    disp('Num of regions in mask does not equal num of clusters to extract from.')
    disp('The most likely cause is an ill-formed/outdated obj.volInfo.cluster field');
    disp('in the region-defining object.  Try obj = reparse_contiguous(obj);');
    disp('Stopping in debugger so you can check/debug.');
    keyboard
end

for i = 1:nregions
    imgvec = maskData == u(i);
    
    regiondat = obj.dat(:, imgvec);
    
    if ~isempty(regiondat)
        
        if any(strcmp(varargin, 'pattern_expression'))
            w = maskvals(imgvec);
            
            if pattern_norm
                % L2 norm of pattern, norm to length 1 -- geometric mean
                %w = w ./ sum(abs(w).^2)^(1/2);
                
                % L1 norm
                w = w ./ sum(abs(w));
                
                cl(i).val_descrip = 'Mask weights, L1 normed';
            else
                cl(i).val_descrip = 'Mask weights, no normalization';
            end
            
            % Calculate similarity, pass in optional inputs for sim metrics
            regionmean = canlab_pattern_similarity(regiondat', w, 'ignore_missing', varargin{:});
            
            %regionmean = regiondat * w;
            
            cl(i).val = w;
            
            if nargout > 1
                % optional outputs:
                % separate pexp for mean across region (unit vector) and
                % mean-centered pattern
                z = ones(size(w));
                w2 = w - mean(w);
        
                clroimean(i).val = z;
                clroimean(i).dat = regiondat * z;
                
                clpattern(i).val = w - mean(w);
                clpattern(i).dat = canlab_pattern_similarity(regiondat', w2, 'ignore_missing', varargin{:});
                
            end
            
        else
            % else, simple average
            
            if size(regiondat, 2) == 1
                regionmean = double(regiondat);
            else
                regionmean = double(nanmean(regiondat,2));
            end
        end % pattern or average...
        
        cl(i).all_data = single(regiondat);
        
    else
        regionmean = NaN .* zeros(size(dat, 1), 1);
    end
    
    cl(i).dat = regionmean;
    
end

if doverbose
    fprintf('Done.\n');
end

% if nargout > 1
%     varargout{1} = clroimean;
%     varargout{2} = clpattern;
% end


end % function
