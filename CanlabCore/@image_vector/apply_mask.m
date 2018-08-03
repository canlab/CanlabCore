function [dat, mask] = apply_mask(dat, mask, varargin)
% Apply a mask image (image filename or fmri_mask_image object) to an image_vector object
% stored in dat.
%
% This can be used to:
%  - Mask an image_vector or fmri_data object with a mask
%  - Obtain "pattern expression" for a weight map (entered as the
%    mask, here) in a series of images stored in dat.
%
% The mask or weight map does not have to be in the same space as the dat;
% it will be resampled to the space of the data in dat.
%
% To extract pattern expression values for each ROI within a mask use extract_roi_averages()
%
% :Optional Inputs:
%
%   **pattern_expression:**
%        calculate and return the dot product of each
%        image in dat and the values in the mask.  This is useful if comparing
%        expression values that are comprised of different datasets or differing
%        number of voxels.
%
%   **correlation:**
%        calculate the pearson correlation coefficient of each
%        image in dat and the values in the mask.
%
%   **norm_mask:**
%        normalize the mask weights by L2 norm, for patt expression
%        only.
%
%   **ignore_missing:**
%        use with pattern expression only. Ignore weights on voxels
%        with zero values in test image. If this is not entered, the function will
%        check for these values and give a warning.
%
%   **invert:**
%        Invert the mask so that out-of-mask voxels are now in (using
%        the mask as an 'exclude mask' rather than an include-mask. If pattern
%        expression is requested, the behavior is different, and it inverts the
%        sign of in-mask pattern weights.
%
%   **cosine_similarity:**
%        use with pattern expression only. scales expression by product of
%        l2 norms (norm(mask)*norm(dat))
%
% :Examples:
% ::
%
%     [dat, mask] = apply_mask(dat, mask)
%     [dat, mask] = apply_mask(dat, mask image name)
%     [dat, mask] = apply_mask(dat, mask image vector object)
%     [pattern_exp_values] = apply_mask(dat, weight map image, 'pattern_expression', 'ignore_missing')
%     [pattern_exp_values] = apply_mask(dat, weight map image, 'pattern_expression', 'ignore_missing','correlation')
%
%
% :See also:
%
% extract_roi_averages, to get individual region averages / local pattern expression
% apply_nps, which does whole-pattern and local regional expression
%
% ..
%    Notes:
%    Last modified: 10/30/11 to add support for masks that are weight maps
%    12/15/13:  Luke Chang - added correlation option for pattern-expression
%    5/10/2016: Phil Kragel - added cosine similarity
%    8/2/2018: Tor Wager - if mask is a statistic_image object, apply .sig field
%
% ..

dopatternexpression = 0; % set options
donorm = 0;
doinvert = 0;

% Now handled by canlab_pattern_similarity
% doignoremissing = 0;
% docorr = 0; %run correlation instead of dot-product for pattern expression
% docosine = 0;

if any(strcmp(varargin, 'pattern_expression'))
    dopatternexpression = 1;
    
    % Now handled by canlab_pattern_similarity
    %     if any(strcmp(varargin, 'ignore_missing'))
    %         doignoremissing = 1;
    %     end
    %
    %     if any(strcmp(varargin, 'cosine_similarity'))
    %         docosine = 1;
    %     end
    %
    %
    %     if any(strcmp(varargin, 'correlation')) % run correlation instead of dot-product
    %         docorr = 1;
    %     end
    
end

if any(strcmp(varargin, 'invert'))
    doinvert = 1;
end

if any(strcmp(varargin, 'norm_mask')) % only good for pattern expression
    donorm = 1;
end

% create mask_image object if we have a filename
if ischar(mask)
    mask = fmri_mask_image(mask);
end

isdiff = compare_space(dat, mask);

if isdiff == 1 || isdiff == 2 % diff space, not just diff voxels
    % == 3 is ok, diff non-empty voxels
    
    % Both work, but resample_space does not require going back to original
    % images on disk.
    %mask = resample_to_image_space(mask, dat);
    mask = resample_space(mask, dat);
    
    % tor added may 1 - removed voxels was not legal otherwise
    %mask.removed_voxels = mask.removed_voxels(mask.volInfo.wh_inmask);
    % resample_space is not *always* returning legal sizes for removed
    % vox? maybe this was updated to be legal
    
    if length(mask.removed_voxels) == mask.volInfo.nvox
        disp('Warning: resample_space returned illegal length for removed voxels. Fixing...');
        mask.removed_voxels = mask.removed_voxels(mask.volInfo.wh_inmask);
    end
    
end

% dat = remove_empty(dat); Tor edited 9/10/17 for computational efficiency

% nonemptydat: Logical index of voxels with valid data, in in-mask space
nonemptydat = get_nonempty_voxels(dat);

dat = replace_empty(dat);

% Check/remove NaNs. This could be done in-object...
mask.dat(isnan(mask.dat)) = 0;

% Replace if necessary
mask = replace_empty(mask);

% if mask is a statistic_image object, apply .sig field
if isa(mask, 'statistic_image')
    
    mask.dat(~mask.sig) = 0;  % for masking purposes
    
end

if doinvert
    
    if dopatternexpression
        mask.dat = -mask.dat;
        disp('Reversing sign of mask values.');
        
    else
        disp('Inverting mask: Out-of-mask areas are now in.');
        maskdat = reconstruct_image(mask);
        maskdat = double(~maskdat);
        mask = rebuild_volinfo_from_dat(mask, maskdat(:));
        mask = resample_space(mask, dat);
    end
    
end

% save which are in mask, but do not replace with logical, because mask may
% have weights we want to preserve
inmaskdat = logical(mask.dat);


% Remove out-of-mask voxels
% ---------------------------------------------------

% mask.dat has full list of voxels
% need to find vox in both mask and original data mask

if size(mask.volInfo.image_indx, 1) == size(dat.volInfo.image_indx, 1)
    n = size(mask.volInfo.image_indx, 1);
    
    if size(nonemptydat, 1) ~= n % should be all vox OR non-empty vox
        nonemptydat = zeroinsert(~dat.volInfo.image_indx, nonemptydat);
    end
    
    if size(inmaskdat, 1) ~= n
        inmaskdat = zeroinsert(~mask.volInfo.image_indx, inmaskdat);
    end
    
    inboth = inmaskdat & nonemptydat;
    
    % List in space of in-mask voxels in dat object.
    % Remove these from the dat object
    to_remove = ~inboth(dat.volInfo.wh_inmask);
    
    to_remove_mask = ~inboth(mask.volInfo.wh_inmask);
    
elseif size(mask.dat, 1) == size(dat.volInfo.image_indx, 1)
    
    % mask vox are same as total image vox
    nonemptydat = zeroinsert(~dat.volInfo.image_indx, nonemptydat);
    inboth = inmaskdat & dat.volInfo.image_indx & nonemptydat;
    
    % List in space of in-mask voxels in dat object.
    to_remove = ~inboth(dat.volInfo.wh_inmask);
    
    to_remove_mask = ~inboth(mask.volInfo.wh_inmask);
    
elseif size(mask.dat, 1) == size(dat.volInfo.wh_inmask, 1)
    % mask vox are same as in-mask voxels in dat
    inboth = inmaskdat & dat.volInfo.image_indx(dat.volInfo.wh_inmask) & nonemptydat;
    
    % List in space of in-mask voxels in .dat field.
    to_remove = ~inboth;
    
    to_remove_mask = ~inboth;
    
else
    fprintf('Sizes do not match!  Likely bug in resample_to_image_space.\n')
    fprintf('Vox in mask: %3.0f\n', size(mask.dat, 1))
    fprintf('Vox in dat - image volume: %3.0f\n', size(dat.volInfo.image_indx, 1));
    fprintf('Vox in dat - image in-mask area: %3.0f\n', size(dat.volInfo.wh_inmask, 1));
    disp('Stopping to debug');
    keyboard
end

dat = remove_empty(dat, to_remove);
mask = remove_empty(mask, to_remove_mask);

if dopatternexpression
    %mask = replace_empty(mask); % need for weights to match
    
    weights = double(mask.dat); % force double b/c of matlab instabilities
    dat.dat = double(dat.dat); % force double b/c of matlab instabilities
    
    if donorm, weights = weights ./ norm(weights); end
    
    % Pass similarity metric in to canlab_pattern_similarity
    similarity_output = canlab_pattern_similarity(dat.dat, weights, varargin{:});
    
    dat = similarity_output;
    
    % CHECK and weight
    % ---------------------------------------------------------
    
    % Pattern weights can have zeros, which may be valid values in voxels,
    % i.e., with binary masks
    % Images with values of zero or NaN are considered out-of-mask, as they
    % are not valid values. These should be excluded from both image and
    % mask when calculating similarity.
    % Thus, there is an asymmetry between pattern mask and image data
    % in considering which voxels to use.
    % Otherwise, all dot product and similarity metrics are standard.
    
    % We also need to calculate bad values on an image-by-image basis, not
    % relying on remove_empty to exclude voxels with ineligible values
    % across the entire set of images.
  
        
    % End check and weight ---------------------------------------------------------
    
    
end % Pattern expression


% Handle special types
% So far, atlas

if isa(dat, 'atlas')
   
    dat = check_properties(dat, 'compress_index');
    
end

end  % Main function




function nonemptydat = get_nonempty_voxels(dat)
empty_voxels = all(dat.dat' == 0 | isnan(dat.dat'), 1)';

if size(empty_voxels, 1) == dat.volInfo.n_inmask
    % .dat is full in-mask length, we have not called remove_voxels or there are none to remove
    nonemptydat = ~empty_voxels;
    
elseif length(dat.removed_voxels) == dat.volInfo.n_inmask
    % .dat is not in-mask length, and we have .removed_voxels defined
    nonemptydat = false(dat.volInfo.n_inmask, 1);
    nonemptydat(~dat.removed_voxels) = true;
    
    % additional: we could have invalid voxels that have been changed/added since
    % remove_empty was last called.
    %     dat.removed_voxels(dat.removed_voxels) =
    %     nonemptydat(empty_voxels
end

end

