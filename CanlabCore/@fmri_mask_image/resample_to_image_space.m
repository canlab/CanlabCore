function obj = resample_to_image_space(obj, sampleto, varargin)
% Resamples data in an fmri_mask_image object (obj) to the space of another
% image (e.g., a functional image, for data extraction)
% The volInfo field will be the same as the sampleto volume info.
% The mask will have zeros in obj.dat for out-of-mask voxels.
% THIS FUNCTION USES SCN_MAP_IMAGE AND REQUIRES THAT THE ORIGINAL IMAGE BE
% AVAILABLE ON DISK.  Multiple resamplings will break the function because
% the new space will be different from the original one on disk.  Use the
% more general resample_space.
%
% NOTE: Mask is *reloaded* from original data if space is remapped, and you
% cannot use manual thresholding of the mask. This is a feature of the
% map_to_image_space method and scn_map_image
%
% :Usage:
% ::
%
%     obj = resample_to_image_space(obj, sampleto <img name or image_vector object>)
%
% :Inputs:
%
%   **obj:**
%        must be an fmri_mask_image object
%   **sampleto:**
%        can be either:
%           1. An image name to sample to
%           2. Another fmri_mask_image object (but image must exist on path!)
%
% :Optional inputs:
%   **mask:**
%        Apply sampleto as mask so that only voxels in the sampleto mask
%        are retained in obj.dat.
%
% THIS FUNCTION WORKS, BUT IS DEPRECATED BECAUSE RESAMPLE_SPACE IS MORE
% GENERAL.  resample_space does not require the resampling of the original
% image from disk, which this does.  resample_space is slower, though.
% 
% See Also: resample_space, for a method that does not require images to
% exist on disk on the path.
%
% ..
%    Programmers' notes
%    Tor: July 2011: Edited because old version will apply mask when resampling
%    space. Edited default behavior to NOT mask with voxels only in sampleto
%    space.  This could cause bugs in other functions that need to be worked
%    out.  The optional argument 'mask' should produce the old default
%    behavior.
%
%    Oct 30, 2011: obj.dat field after sampling did not conform to standard,
%    because only in-mask voxels in volInfo were not selected.  This was
%    fixed.
% ..

switch class(sampleto)
    case 'char'
        image_name_to_sample_to = sampleto;
        
        volInfo_to = iimg_read_img(sampleto, 2, 1, 1); % read data from file, first volume only
        
    case {'image_vector', 'fmri_data', 'fmri_mask_image', 'statistic_image'}
        
        if any(strmatch('space_defining_image_name', fieldnames(sampleto), 'exact')) && ~isempty(sampleto.space_defining_image_name)
            image_name_to_sample_to = sampleto.space_defining_image_name;
        else
            image_name_to_sample_to = sampleto.volInfo.fname;
        end
        
        volInfo_to = sampleto.volInfo;
        
    otherwise
        error('fmri_mask_image.resample_to_image_space: illegal sampleto input.');
end

obj.space_defining_image_name = image_name_to_sample_to;

% do most of the work here.
obj.dat = scn_map_image(obj, volInfo_to);

obj.dat = obj.dat(:);

obj.history{end+1} = sprintf('Used image %s', obj.volInfo.fname);
obj.history{end+1} = sprintf('Sampled to space of %s', image_name_to_sample_to);



obj.dat(isnan(obj.dat)) = 0;

if any(strcmp(varargin, 'mask'))
    obj.volInfo = volInfo_to; % This will define mask with voxels in sampleto image

    % Need the same voxels in each
    % Possibilities: mask could have voxels outside image-defining mask
    if size(obj.dat, 1) == volInfo_to.n_inmask % ~isfield(volInfo_to, 'n_inmask') ||
        % OK, do nothing
    elseif size(obj.dat, 1) == size(volInfo_to.image_indx, 1)
        obj.dat = obj.dat(volInfo_to.image_indx, :);
    else
        disp('Error: mask data is wrong size. Debug me.')
        keyboard
    end
else

    % Update volInfo structure to reflect resampling
    % -------------------------------------------------------

    obj.volInfo.fname = 'REMOVED: CHANGED SPACE';
    obj.volInfo.mat = volInfo_to.mat;
    obj.volInfo.dim = volInfo_to.dim;
    obj.volInfo.dt = volInfo_to.dt;
    
    
    obj.volInfo(1).descrip = sprintf('Space of %s', image_name_to_sample_to);
    
    % Stuff for extended volInfo entries
    obj.volInfo(1).image_indx = ~isnan(obj.dat) & abs(obj.dat) > 10*eps ;    % indices of all voxels in mask
    obj.volInfo(1).nvox = length(obj.dat);
    
    obj.volInfo(1).wh_inmask = find(obj.volInfo(1).image_indx);    % indices of all voxels in mask
    obj.volInfo(1).n_inmask = length(obj.volInfo(1).wh_inmask);
    
    [i, j, k] = ind2sub(obj.volInfo(1).dim(1:3), obj.volInfo(1).wh_inmask);
    obj.volInfo(1).xyzlist = [i j k];
    
    if obj.volInfo(1).n_inmask < 50000
        obj.volInfo(1).cluster = spm_clusters(obj.volInfo(1).xyzlist')';
    else
        obj.volInfo(1).cluster = ones(obj.volInfo(1).n_inmask, 1);
    end

    % obj.dat should contain only in-mask values.
    obj.dat = obj.dat(obj.volInfo(1).wh_inmask, :);
    
end



obj.removed_voxels = false(size(obj.dat, 1), 1);


end % function

