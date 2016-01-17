function voldata = iimg_reconstruct_vols(dat, volInfo, varargin)
% Reconstruct a 3-D or 4-D volume from a "dat" matrix of vectorized images
%
% :Usage:
% ::
%
%     voldata = iimg_reconstruct_vols(dat, volInfo, [optional args])
%
% :Optional Inputs: If entered, will write .img file to disk
%
%   **'outname':**
%        followed by output name
%
%   **'descrip':**
%        followed by description for img file
%
%   **'slice':**
%        followed by slice number of single-slice data in image
%
% :Example of image reading:
%
%   1) Get volume info from first volume of a 3-D or 4-D image
%       volInfo structure has necessary info for converting to/from
%       "vectorized" format
%       dat returns 4-D data from entire image (all volumes)
%       ::
%
%           img_name = 'test_run1_pca.img';
%           [maskInfo, dat] = iimg_read_img(img_name, 2);
%
%   2) Get data in "vectorized" image format for each volume in the
%      image. Works for a list of images too. data is in-mask voxels x volumes
%      this can be useful when you want to return whole-brain data for many
%      images, but in a search volume only (i.e., no extra-brain voxels)
%      ::
%
%            data = iimg_get_data(maskInfo, img_name);
% data = data';  % make sure columns are volumes
%
%   3) Write out a 4-D image with the same data, called test_run1_pca2.img
%      ::
%
%           voldat3D = iimg_reconstruct_vols(data, maskInfo, 'outname', 'test_run1_pca2.img');
%
% ..
%    NOTES:
%    2013/3/25: 
%    Luke[ea] added optional input to retain original datatype (default = float32)
% ..

    is_writing_file = 0;
    is_single_slice = 0;
    keepdt = 0;
    slice_number = NaN;
    descrip = 'Created by iimg_reconstruct_vols';

    % This stuff is now handled by scn_write_plane...
    % Do not force spm scaling factors and get rid of private, which may
    % contain information from other images/irrelevant info for this...
    % %     if isfield(volInfo, 'pinfo'), volInfo = rmfield(volInfo, 'pinfo'); end
    % %     if isfield(volInfo, 'private'), volInfo = rmfield(volInfo, 'private'); end

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'outname'
                    outname = varargin{i+1};
                    is_writing_file = 1;

                case 'descrip'
                    descrip = varargin{i+1};

                case 'slice'
                    is_single_slice = 1;
                    slice_number = varargin{i+1};
                    
                case 'keepdt'
                    keepdt = 1;
            end
        end
    end

    % strip scaling factors from volInfo so that they will be
    % recalculated in a way appropriate for this data.
    % volInfo is intended to be dim info, but not scaling factor info.

    if isfield(volInfo, 'pinfo'), volInfo = rmfield(volInfo, 'pinfo'); end

    % check data and make sure it conforms to full-length index vector
    % adjust if only in-mask voxels are in data
    % added dtfloat as optional input. - added 3/25/13 luk[ea]
    [dat volInfo] = check_dat(dat, volInfo, is_single_slice, keepdt, slice_number );

    nimgs = size(dat, 2); % one column in vectorized format per 3D image

    % Create empty array
    if is_single_slice
        voldata = zeros([volInfo.dim(1:2) nimgs]);
    else
        voldata = zeros([volInfo.dim(1:3) nimgs]);
    end

    % Fill with data
    for i = 1:nimgs
        if is_single_slice
            voldata(:, :, i) = reshape(dat(:, i), volInfo.dim(1:2));
        else
            voldata(:, :, :, i) = reshape(dat(:, i), volInfo.dim(1:3));
        end
    end

    if is_writing_file
        if is_single_slice
            write_file(voldata, volInfo, outname, descrip, slice_number);
        else
            write_file(voldata, volInfo, outname, descrip);
        end
    end

end  % End Main Function


% ------------------------------------------------------------------
% * Check and make sure data is in full list; adjust if mask-only; alter datatype of volInfo if needed
% ------------------------------------------------------------------
function [dat volInfo] = check_dat(dat, volInfo, is_single_slice, keepdt, varargin)
    % last arg is slice no., if applicable

    if ~isfield(volInfo, 'nvox')
        error('volInfo.nvox and other fields not found.  Use iimg_read_image.m to prepare volInfo.');
    end

    ndat = size(dat, 1);

    if is_single_slice
        % data should either equal prod(volInfo(1:2)) or fewer elements, in
        % which case one slice will be reconstructed
        nslice = prod(volInfo.dim(1:2));
        %z = volInfo.xyzlist(:,3);
        slice_number = varargin{1};

        if ndat ~= nslice
            % get indices in full volume of this slice
            st = 1 + (slice_number - 1) .* nslice;  % starting index of this slice
            en = slice_number .* nslice;

            tmpdat = zeros(nslice, 1); % full volume, so we know how voxels map into 3-D array
            tmpdat(volInfo.image_indx(st:en)) = dat;
            dat = tmpdat;
        end


    elseif ndat ~= volInfo.nvox
        % if data already includes all voxels, do nothing, otherwise assume data is only in-mask, and reconstruct
        if ~isfield(volInfo, 'n_inmask')
            error('volInfo.n_inmask not found.  Use iimg_read_img.m with extended output flag to prepare volInfo.');
        end

        if ndat ~= volInfo.n_inmask
            error('volume info in struct does not seem to match number of voxels in index image.');
        else
            % we have in-mask data only and have to put it back in the right
            % place in the original image
            nimgs = size(dat, 2); % number of images
            tmpdat = zeros(volInfo.nvox, nimgs);
            for i = 1:nimgs
                tmpdat(volInfo.image_indx, i) = dat(:, i);
            end
            dat = tmpdat;
        end
    end

    % alter data type if needed
    switch(spm('Ver'))
        case 'SPM2'
            if spm_type(volInfo.dim(4), 'swapped')
                disp('Swapped datatypes are supported by SPM2 but not SPM5. Un-swapping.');
                volInfo.dim(4) = volInfo.dim(4) ./ 256;

                if(spm_type(volInfo.dim(4), 'intt'))
                    volInfo.dim(4) = spm_type('float');
                end
            end

        case {'SPM5', 'SPM8', 'SPM12'} %added keepdt to retain original dt if requested : luk(ea)
            %Outputs float32 by default otherwise keeps original data type
            %if 'keepdt' is used as optional input.
            if isfield(volInfo, 'dt') && spm_type(volInfo.dt(1), 'intt') && ~keepdt
               volInfo.dt(1) = spm_type('float32');
            end
            
        otherwise 
            error('Unrecognized SPM type.  Please update code or use SPM2/5/8/12!');
    end
end

% ------------------------------------------------------------------
% Write full image or slice, if last arg. entered
% ------------------------------------------------------------------
function write_file(voldata, volInfo, outname, descrip, varargin)

    sz = size(voldata);
    % get number of images for volume or slice data
    if length(sz) > 3
        nimgs = sz(4);
    elseif ~isempty(varargin)
        % we have one slice for nimgs images, slices are concat. in 3rd dim
        if length(sz) == 2  % we've passed in only a slice
            nimgs = 1;
        else
            nimgs = sz(3);
        end
    else
        nimgs = 1;
    end

    if ~isempty(varargin)
        % We have a slice, and use spm_write_plane.m
        % write one slice for a series of images
        slice_number = varargin{1};

        filenames = [];
        for i = 1:nimgs
            filenames = char(filenames, make_img_filename(outname, i));
        end
        filenames = filenames(2:end, :);

        % Note: pinfo is set to [1 0 0] for new images; no rescaling of
        % images.
        scn_write_plane(filenames, voldata, slice_number, volInfo);

    else
        warning('off'); % empty images return many warnings
        for i = 1:nimgs
            [volInfo.fname, volInfo.n(1)] = make_img_filename(outname, i);
            volInfo.descrip = descrip;

            % Note: done above: strip scaling factors from volInfo so that they will be
            % recalculated in a way appropriate for this data.
            % volInfo is intended to be dim info, but not scaling factor info.

            spm_write_vol(volInfo, squeeze(voldata(:, :, :, i)));
        end
        warning('on');
    end
end


function [name, n] = make_img_filename(name, imagenum)
    [d, name, ext] = fileparts(name);
    name(name == '.') = '_';

    % If you pass in a filename with a ,### suffix, it will keep that suffix by
    % using the volInfo.n field, which tells spm_write_vol which image to write
    % to
    % If you pass in a single name and are writing 4-D images, it will
    % increment the image number counter.
    %
    % changed nov 07 to preserve commas for multi-volume images.
    % handle the comma to index volume: return as n and take it off the
    % filename
    t = find(ext==',');

    if ~isempty(t)
        t = t(1);
        n1 = ext((t+1):end);
        if ~isempty(n1),
            n = str2num(n1);
            ext = ext(1:(t-1));
        end
    else
        n = imagenum;
    end

    name = fullfile(d, [name ext]);
end
