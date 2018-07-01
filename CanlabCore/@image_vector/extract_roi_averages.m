function cl = extract_roi_averages(obj, mask, varargin)
% This image_vector method a extracts and averages data stored in an fmri_data object
% from a set of ROIs defined in a mask.
% It is *slightly* different from the fmri_data method, as fmri_data has
% more fields.
%
% This version resamples the mask_image to be in the same space as the obj if needed.
%
% Regions to average over can be either regions of contiguous voxels
% bounded by voxels with values of 0 or NaN, which are considered non-data
% values, or regions defined by unique integer codes in the mask image
% (i.e., for atlas images with unique codes for each defined region.)
%
% Mask/Atlas image does NOT have to be in the same space as the images to
% extract from.  It will be remapped/resliced.
%
% Extracted data is returned in single data format.
%
% :Usage:
% ::
%
%    cl = extract_roi_averages(image_vector obj, mask, [average_over])
%
% :Inputs:
%
%   - char array of strings containing 4D image file names (data extracted from these)
%   - mask_image to extract from.
%
% :Optional inputs:
%
%   **average_over:**
%      - Default: 'contiguous_regions' to average over contiguous voxels
%        bounded by voxels of 0 or NaN (non-data values)
%      - Alt. option = 'unique_mask_values' to average over unique integer codes in the mask image
%        (i.e., for atlas images with unique codes for each defined region)
%
% :Examples:
% ::
%
%    imgs_to_extract_from = filenames('w*.nii','char');
%    mask_image = which('anat_lbpa_thal.img');
%    [cl, imgdat] = extract_image_data(imgs_to_extract_from, mask_image);
%
% :See also:
%
% For an non-object-oriented alternative, see extract_image_data.m
% fmri_data.extract_roi_averages, apply_parcellation.m

% Programmers' notes:
% This function is different from fmri_data.extract_roi_averages
% Better to have only one function of record in the future...
% Note: 
% cl = extract_roi_averages(imgs, atlas_obj{1});
% accomplishes the same task as apply_parcellation, returns slightly different values due to interpolation.  


doverbose = 1;
if any(strcmp(varargin, 'noverbose'))
    doverbose = 0;
end

if doverbose
    fprintf('image_vector.extract_roi_averages: ');
end

if ~isa(mask, 'image_vector')
    mask = image_vector('image_names', mask);
end

obj = replace_empty(obj);
mask = replace_empty(mask);

% --------------------------------------------
% Resample space if needed
% --------------------------------------------

isdiff = compare_space(obj, mask);

if isdiff == 0 || isdiff == 3 % Same space, diff voxels
   % do nothing
elseif isdiff == 2
    error('Invalid object: volInfo structure missing for one or more objects');
else
    mask = resample_space(mask, obj);
end

% --------------------------------------------
% Redefine elements of new mask to apply
% to keep only coords in original dataset
% --------------------------------------------

is_inmask = obj.volInfo.image_indx & mask.volInfo.image_indx;
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
wh_to_keep = is_inmask(obj.volInfo.wh_inmask);

% eliminate out-of-mask voxels before indexing into them with new mask
obj.dat = obj.dat(wh_to_keep, :)';

if doverbose
    fprintf('\n');
end

% ---------------------------------
% define region object based on choices
% ---------------------------------

average_over = 'unique_mask_values'; %'contiguous_regions'  or 'unique_mask_values';

for varg = 1:length(varargin)
    if ischar(varargin{varg})
        switch varargin{varg}
            
            % reserved keywords
            case 'contiguous_regions', average_over = 'contiguous_regions';
            case 'unique_mask_values', average_over = 'unique_mask_values';
                
            case 'pattern_expression', error('pattern expression only defined for fmri_data.extract_roi_averages.');
                
            case 'noverbose'  % do nothing; already done
                
            otherwise
                disp('fmri_data.extract_roi_averages: Illegal string value for average_over.');
                fprintf('You entered ''%s''\n Valid values are %s or %s\n', varargin{varg}, '''contiguous_regions''', '''unique_mask_values''');
                error('Exiting');
        end
    end
end

if isa(mask, 'atlas') % If atlas object
    cl = atlas2region(mask);                % handles a few things better
else
    cl = region(mask, average_over);
end

cl(1).source_images = obj.fullpath;

% ---------------------------------
% Now get averages by cluster
% ---------------------------------

if doverbose
    fprintf('Averaging data. ');
end

switch size(mask.dat, 1)
    case length(mask.volInfo.wh_inmask)
        %maskData = mask.dat(logical(mask.volInfo.wh_inmask), :);
        % We already have in-mask data; do not select
        maskData = mask.dat;
        
    case length(mask.volInfo.image_indx)
        maskData = mask.dat(logical(mask.volInfo.image_indx), :);
    otherwise
        error('mask .dat is wrong size.')
end

switch average_over
    
    % Define integer codes for sets of voxels to average over.
    
    case 'unique_mask_values'
        maskData = round(maskData);
        u = unique(maskData)'; u(u == 0) = [];
        nregions = length(u);
        fprintf('Averaging over unique mask values, assuming integer-valued mask: %3.0f regions\n', nregions);
        
        
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
        if size(regiondat, 2) == 1
            regionmean = double(regiondat);
        else
            regionmean = double(nanmean(regiondat')');
        end
        cl(i).all_data = single(regiondat);
        
    else
        regionmean = NaN .* zeros(size(dat, 1), 1);
    end
    
    cl(i).dat = regionmean;
    
end

if doverbose
    fprintf('Done.\n');
end

end % function
