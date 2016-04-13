function [img, Vto] = scn_resample_voxel_size(loadImg, voxsize, varargin)
% This function takes an image name in loadImg
% and loads the data, resampling to the space defined
% in the image sampleTo.
%
% :Usage:
% ::
%
%     [imgData, volInfo_mapto] = scn_resample_voxel_size(loadImg, voxsize, varargin)
%
% take volInfo and fmri_mask_image inputs as well
% as image file names
%
% :Optional Inputs:
%
%   **'write':**
%        followed by name of resampled image to write
%
% Compatible with SPM5/8.
%
% Input images can have the following formats:
%   1) String with name of image file (.img or .nii)
%   2) spm_vol-style V struct (see spm_vol)
%   3) volInfo struct (see iimg_read_img)
%   4) fmri_mask_image object (see fmri_mask_image)
%
% :Examples:
% ::
%
%    % Reslice standard brain mask to 3 x 3 x 3 voxels.
%    img = which('brainmask.nii');
%    [dat, Vto] = scn_resample_voxel_size(img, [3 3 3], 'write', 'test.img');
%    spm_image('init', 'test.img');
%    spm_check_registration(char(img, 'test.img'));
%
% ..
%    Tor Wager, Oct 2010
% ..

    % ..
    %    optional input arguments
    % ..
    if ~isempty(varargin)
        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}

                    % functional commands
                    case {'write', 'name', 'outname'}
                        outname = varargin{i+1};
                        varargin{i+1} = [];

                    otherwise, warning(['Unknown input string option:' varargin{i}]);
                end
            end
        end
    end


    % image whose space to sample into
    % -------------------------------------------------------
    sampleTo = loadImg;  % we're just going to change the voxel size
    
    if ischar(sampleTo)
        sampleTo = deblank(sampleTo);
        Vto = spm_vol(sampleTo);               % volume to sample TO
        
    % Could already be a volInfo structure
    elseif strcmp(class(sampleTo), 'fmri_mask_image')
        Vto = sampleTo.volInfo;
        
    elseif isstruct(sampleTo)
        Vto = sampleTo;
    end
    
    Vto = Vto(1); % for SPM5 compat
    
    if isempty(Vto)
        fprintf('%s is empty or missing.\n', sampleTo)
        return
    end

    Mto = Vto(1).mat;                        % mat to sample TO

    % now apply voxel size change
    for i = 1:3
        Mto(i, i) = sign(Mto(i, i)) * voxsize(i);
    end

    % image to load and sample to new space
    % -------------------------------------------------------
    if ischar(loadImg)
        loadImg = deblank(loadImg);
        Vmap = spm_vol(loadImg);                % mask to resample
    
    % Could already be a volInfo structure
    elseif strcmp(class(loadImg), 'fmri_mask_image') | strcmp(class(loadImg), 'image_vector') | strcmp(class(loadImg), 'fmri_data')
        Vmap = loadImg.volInfo;
        
    elseif isstruct(loadImg)
        Vmap = loadImg;
        
    end
        
    if isempty(Vmap)
        fprintf('%s is empty or missing.\n', loadImg)
        return
    end

    Mmap = Vmap(1).mat;



    % image dimensions
    % -------------------------------------------------------
    %dim = Vto(1).dim(1:3);
    % re-make dim based on new voxel size
    dim = round(abs(Vto.dim .* diag(Vto.mat(1:3, 1:3))') ./ voxsize);

    % output image data
    % -------------------------------------------------------
    img = zeros(dim);


    % map slice-by-slice
    % -------------------------------------------------------
    for j = 1:dim(3)

        Mslice  = spm_matrix([0 0 j]);      % Matrix specifying this slice

        Mtrans  = Mto \ Mmap \ Mslice;          % Affine mappping mtx: Mask -> TOvol

        img(:, :, j) = spm_slice_vol(Vmap(1), Mtrans, dim(1:2), [0 NaN]);


    end

    % create image, if asked for
    % -------------------------------------------------------
    if exist('outname', 'var')
        Vout = Vto;
        Vout.fname = outname;
        Vout.mat = Mto;
        Vout.dim = dim;
        Vout  = spm_create_vol(Vout);

        spm_write_vol(Vout, img);
    end



end
