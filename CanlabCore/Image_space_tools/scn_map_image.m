function [img, Vto] = scn_map_image(loadImg, sampleTo, varargin)
% This function takes an image name in loadImg
% and loads the data, resampling to the space defined
% in the image sampleTo. The resampled image will retain
% the data type of the input image.
%
% :Usage:
% ::
%
%     [imgData, volInfo_mapto] = scn_map_image(loadImg, sampleTo, varargin)
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
%    img = scn_map_image(EXPT.mask, EXPT.SNPM.P{1}(1,:), 'write', 'resliced_mask.img');
%
% ..
%    Tor Wager, July 2007
%    edited, Oct 2010, to take volInfo and fmri_mask_image inputs as well
%    as image file names
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
    if ischar(sampleTo)
        sampleTo = deblank(sampleTo);
        Vto = spm_vol(sampleTo);               % volume to sample TO
        
    % Could already be a volInfo structure
    elseif isa(sampleTo, 'image_vector')
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


    % image to load and sample to new space
    % -------------------------------------------------------
    if ischar(loadImg)
        loadImg = deblank(loadImg);
        Vmap = spm_vol(loadImg);                % mask to resample
    
    % Could already be a volInfo structure
    elseif isa(loadImg, 'image_vector') % also includes 'fmri_mask_image' 'fmri_data'
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
    dim = Vto(1).dim(1:3);


    % output image data
    % -------------------------------------------------------
    img = zeros(dim);


    % check for space change
    % -------------------------------------------------------
    if strcmp(Vmap(1).fname, 'REMOVED: CHANGED SPACE')
        
        disp('You have re-mapped the space of this image, and cannot use scn_map_image.m')
        disp('A common reason for this error is trying to map a data object to another space twice,')
        disp('which cannot currently be done.');
        error('Unsupported operation: See message above.')
        
    end
    
    
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
        Vout.dt(1) = Vmap.dt(1); % keep data type of input image
        Vout.fname = outname;
        Vout  = spm_create_vol(Vout);

        spm_write_vol(Vout, img);
    end



end
