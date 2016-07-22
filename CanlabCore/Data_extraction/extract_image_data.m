function [imgdat, volInfo, cl] = extract_image_data(imgs_to_extract_from, mask_image, varargin)
% Generic function for extracting image data from a mask or atlas image,
% and returning the data and averages within regions specified by the user.
%
% :Usage:
% ::
%
%    [imgdat, volInfo, cl] = extract_image_data(imgs_to_extract_from, mask_image, varargin)
%
% Regions to average over can be either regions of contiguous voxels
% bounded by voxels with values of 0 or NaN, which are considered non-data
% values, or regions defined by unique integer codes in the mask image
% (i.e., for atlas images with unique codes for each defined region.)
% 
% Mask/Atlas image does NOT have to be in the same space as the images to
% extract from.  It will be remapped/resliced.
%
% extracted data is returned in single data format.
%
% :Inputs:
%
%   **char array** of strings containing 4D image file names (data extracted from these)
%
%   **mask_image** to extract from
% 
% :Optional inputs:
%
%   **average_over:**
%        Default = 'contiguous_regions' to average over contiguous voxels
%        bounded by voxels of 0 or NaN (non-data values)
%
%        Alt. option = 'unique_mask_values' to average over unique integer codes in the mask image
%        (i.e., for atlas images with unique codes for each defined
%        region)
% 
% :Examples:
% ::
%
%    imgs_to_extract_from = filenames('w*.nii','char');
%    mask_image = which('anat_lbpa_thal.img');
%    [imgdat, volInfo, cl] = extract_image_data(imgs_to_extract_from, mask_image, 'unique_mask_values');
% 
% :Related functions:
% For an object-oriented alternative, see the fmri_data class and extract_roi_averages method

    
    if nargin < 2 || isempty(mask_image)
        mask_image = which('anat_lbpa_thal.img');
        fprintf('Using default mask: %s\n', mask_image);
        if isempty(mask_image), error('Cannot find mask image!'); end
    end

    average_over = 'contiguous_regions'; %'contiguous_regions'  or 'unique_mask_values';

    for varg = 1:length(varargin)
        if ischar(varargin{varg})
            switch varargin{varg}

                % reserved keywords
                case 'contiguous_regions', average_over = 'contiguous_regions';
                case 'unique_mask_values', average_over = 'unique_mask_values';
            end
        end
    end

    space_defining_image = deblank(imgs_to_extract_from(1, :));


    % Note: we need to write out the image here only because we need its
    % anatomical boundaries (in-mask areas), contig voxels, etc. for later

    tmpname = num2str(round(10*rand(1, 5))); tmpname(tmpname == ' ') = []; 
    tmpname = ['tmp_mask_' tmpname '.img'];


    maskData = scn_map_image(mask_image, space_defining_image, 'write', tmpname);
    maskData = maskData(:);


    volInfo = iimg_read_img(tmpname, 2);
    maskData = maskData(volInfo.wh_inmask);


    % This needs to be made Windows-compatible. kludgy now.
    delete(tmpname);
    delete([tmpname(1:end-4) '.hdr']);

    volInfo.fname = sprintf('%s (resliced to image space)', mask_image);
    volInfo.maskname = mask_image;
    
    % Now we have the mask data and volInfo structure, and we can extract




    %%
    clear imgdat
    switch spm('Ver')


        case {'SPM12', 'SPM8', 'SPM5'}


            imgdat = iimg_get_data(volInfo, imgs_to_extract_from, 'single', 'verbose', 'noexpand');


        case {'SPM2', 'SPM99'}
            % legacy, for old SPM


            imgdat = iimg_get_data(volInfo, imgs_to_extract_from, 'single', 'verbose');


        otherwise
            error('Unknown version of SPM! Update code, check path, etc.');
    end


    %% Now get averages by cluster, if requested


    switch average_over


        % Define integer codes for sets of voxels to average over.


        case 'unique_mask_values'
            maskData = round(maskData);
            u = unique(maskData)'; u(u == 0) = [];
            nregions = length(u);
            fprintf('Averaging over unique mask values, assuming integer-valued mask: %3.0f regions\n', nregions);


        case 'contiguous_regions'
            u = unique(volInfo.cluster); u(u == 0) = [];
            maskData = volInfo.cluster;


        case 'none'
            cl = [];
            return
           

        otherwise
            error('Illegal value for average_over.  See help for this function.');
    end


    %% Now get the average activity in each region and define a "cluster"


    nregions = length(u);


    cl = struct('title', 'title', 'threshold', NaN, 'Z', NaN, 'voxSize', abs(diag(volInfo.mat(1:3, 1:3)))', ...
        'XYZ', [0 0 0]', 'XYZmm', [0 0 0]', 'from_label', NaN, 'M', volInfo.mat, 'dim', volInfo.dim);


    cl(1:nregions) = cl;




    for i = 1:nregions
        imgvec = maskData == u(i);


        regiondat = imgdat(:, imgvec);


        if ~isempty(regiondat)
            regionmean = double(nanmean(regiondat')');


        else
            regionmean = NaN .* zeros(size(dat, 1), 1);
        end

        cl(i).average_data = regionmean;
        cl(i).XYZ = volInfo.xyzlist(imgvec,:)';
        cl(i).XYZmm = voxel2mm(cl(i).XYZ,volInfo.mat);

        cl(i).Z = ones(1,size(cl(i).XYZ,2));

    end
end
