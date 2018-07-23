function [volInfo, dat] = iimg_read_img(inputimgs, extended_output_flag, reading_data, first_vol)
% :Usage:
% ::
%
%     [volInfo, dat] = iimg_read_img(inputimgs, [extended_output_flag (1|2)], [reading_data (0|1)], [first vol only (0|1)]);
%
% - reads 3-D or 4-D img/nii files, with extended volInfo output
% - flag to read first volume only
% - can read iimg_data vectorized form as well as file names (just passes dat to output)
% - uses SPM
%
% volInfo fields:
%   - all fields returned by spm_vol, plus:
%   - .nvox - number of voxels in 3d image
%  - .image_indx - logical index (e.g., [0 1 1 0 1 ...]) of all nonzero voxels
%  - .wh_inmask (extended_output_flag > 0) - linear index (e.g., [2 3 5 ...]) - same as image_indx, but as a result of find
%  - .n_inmask (extended_output_flag > 0) - length of .wh_inmask
%  - .xyzlist (extended_output_flag > 0) - voxel coords of .wh_inmask voxels
%  - .cluster (extended_output_flag > 1) - cluster structure of in-mask voxels
%
% Extended output flag values:
%   1: Add xyz coord list and linear index of which voxels nonzero & non-nan in mask
%   2: Add clusters from spm_clusters
%
% Enforces that dat is a vector of data values, regardless of input.
%
% :Examples:
% ::
%
%    imname = 'rob_tmap_0001_filt_t_3-05_k5_pos.img';
%    volInfo = iimg_read_img(imname);
%
% WARNING: 
% The behavior of this function is SPM version dependent.
% volInfo only contains the spm_vol structure of the FIRST image!!!
% Data in SPM5/8 should be 4-D, but in SPM2 should be 3-D
%
% Example of image reading/writing
% -------------------------------------------------------------------------
%   1) Get volume info from first volume of a 3-D or 4-D image
%       volInfo structure has necessary info for converting to/from
%       "vectorized" format
%       dat returns 4-D data from entire image (all volumes)
%       ::
%
%           iimg_name = 'test_run1_pca.img';
% [maskInfo, dat] = iimg_read_img(img_name, 2);
%
%   2) Get data in "vectorized" image format for each volume in the
%      image. Works for a list of images too. data is in-mask voxels x volumes
%      this can be useful when you want to return whole-brain data for many
%      images, but in a search volume only (i.e., no extra-brain voxels)
%      ::
%
%           idata = iimg_get_data(maskInfo, img_name);
%           data = data';  % make sure columns are volumes
%
% % 3) Write out a 4-D image with the same data, called test_run1_pca2.img
%       ::
%
%           ivoldat3D = iimg_reconstruct_vols(data, maskInfo, 'outname',
%                                                 'test_run1_pca2.img');


% ..
%    Set up arguments
% ..
if ~exist('inputimgs', 'var') || isempty(inputimgs)
    
    error('Input image not found. Check file and path names.');
    
end

if ~exist('extended_output_flag', 'var') || isempty(extended_output_flag)
    extended_output_flag = 0;
end
if ~exist('reading_data', 'var') || isempty(reading_data)
    reading_data = 1;
end

if ~exist('first_vol', 'var') || isempty(first_vol)
    first_vol = 0;
end

if(~reading_data && extended_output_flag)
    warning('If not reading in file data, the extended_output_flag has no meaning');
end

imgtype = get_img_type(inputimgs);

volInfo = [];
num_imgs = size(inputimgs, 1);
if first_vol, num_imgs = 1; end

switch imgtype
    case 'name'
        % ---------------------------------------------------------------
        % * String matrix of filenames: Load first volume, then all data
        % ---------------------------------------------------------------
        % use iimg_read_from_file
        
        space_defining_image = inputimgs(1, :);
        
%         if first_vol
%             whcomma = find(space_defining_image == ',');
%             if ~isempty(whcomma)
%                 space_defining_image = space_defining_image(1:whcomma-1);
%             end
%             
%             space_defining_image = [space_defining_image ',1'];
%             
%         end
     
        % Read first image, or first set of volumes in 4-d file:
        % If firstvol, the ",1" will be added to space defining image after
        % (possible) unzipping and checks on whether it exists.
        
        [volInfo, dat] = iimg_read_from_file(space_defining_image, extended_output_flag, reading_data, first_vol);

        if(reading_data && num_imgs > 1)
            % load other images; check to make sure they're in same dims as first
            
            dat = [dat zeros(volInfo.nvox, num_imgs-1)];

            for i = 2:num_imgs
                [indx data] = iimg_read_from_file(inputimgs(i, :), 0, reading_data, 0); % read all, not just first vol

                if length(data) ~= volInfo.nvox
                    fprintf('Warning: img dims for %s do not match 1st image. Resampling.\n', inputimgs(i, :));
                    
                    data = scn_map_image(inputimgs(i, :), inputimgs(1, :));
                    data = data(:);
                
                end

                dat(:, i) = data;
            end
        end

    case '3dvols'
        % --------------------------------------
        % * matrix of 3-D image data volumes
        % --------------------------------------

        [x, y, z, n] = size(inputimgs);
        nvox = prod(x, y, z);
        dat = zeros(nvox, n);

        for i = 1:num_imgs
            tmp = inputimgs(:, :, :, i);
            dat(:, i) = tmp(:);
        end

        dat(isnan(dat)) = 0;

    case 'indx'
        % --------------------------------------
        % * Already in index -- target format
        % --------------------------------------
        dat = inputimgs;
end
end





function imgtype = get_img_type(inputimgs)
if ischar(inputimgs)
    imgtype = 'name';

    return
end

[x, y, z] = size(inputimgs);

if isa(inputimgs, 'statistic_image')
    imgtype = 'statistic_image';
elseif z > 1
    % 3-D image(s)
    imgtype = '3dvols';
elseif x > 1
    % voxels x images 2-D index matrix
    imgtype = 'indx';
else
    error('I don''t recognize the data format of inputimgs.');
end
end


% iimg_read_from_file
% -------------------------------------------------------------------------

function [volInfo, dat] = iimg_read_from_file(imname, extended_output_flag, reading_data, dofirstvolonly)
% Read image data from file using spm_vol and spm_read_vols.  
% Call vol_file_check, which checks zipped and unzipped versions and
% gunzips if needed

% Returns filename for unzipped image with formatting ,# to read first vol
% only if requested
imname = vol_file_check(imname, dofirstvolonly);

% %  tor: oct 2010: speed up by adding ,1 to get only first vol in spm
% [pth, nam, ext] = fileparts(imname);
% wh_comma = find(ext == ',');
% if isempty(wh_comma)
%     imname = [imname ',1'];
% end

volInfo = spm_vol(imname);

% Changed Feb 2008 to deal with 4-D image reading in SPM5 and above
nimgs = length(volInfo);

for i = 1:nimgs
    volInfo(i).nvox = prod(volInfo(i).dim(1:3));
end

dat = [];

if(reading_data)
    
    switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

        case {'SPM5', 'SPM8', 'SPM12'}
            % spm_defaults is a function
            spm_defaults()

        otherwise
            % unknown SPM
            disp('Unknown version of SPM!');
            spm_defaults()
    end


    % pre-allocate array
    dat = zeros(volInfo(1).nvox, nimgs);

    for i = 1:nimgs
        
        dat_vol = spm_read_vols(volInfo(i));
        dat(:, i) = dat_vol(:);

        volInfo(i).image_indx = dat(:) ~= 0; % locate all voxels in-mask
        
    end

    if extended_output_flag
        % Return this for first image only right now.
        volInfo(1).wh_inmask = find(volInfo(1).image_indx);    % indices of all voxels in mask

        volInfo(1).n_inmask = size(volInfo(1).wh_inmask, 1);

        [i, j, k] = ind2sub(volInfo(1).dim(1:3), volInfo(1).wh_inmask);
        volInfo(1).xyzlist = [i j k];

        if extended_output_flag > 1
            if volInfo(1).n_inmask < 50000
                volInfo(1).cluster = spm_clusters(volInfo(1).xyzlist')';
            else
                volInfo(1).cluster = ones(volInfo(1).n_inmask, 1);
            end
        end
    end

    dat(isnan(dat)) = 0;
end


% Return only the first volInfo struture because it contains generic info
% about all volumes.  Don't need to return every volume for our
% purposes.
volInfo = volInfo(1);

end



function [imname, waszipped] = vol_file_check(imname, varargin)
% Check that filename exists, in either gz (zipped) or unzipped format, and
% return file name of file with full path name. 
% If both zipped and unzipped exist, prefer unzipped.
% Return flag waszipped so file can be re-zipped after loading.

waszipped = false;

if length(varargin) > 0 
    dofirstvolonly = varargin{1};
end

if(~ischar(imname) || isempty(imname))
    error('Image name is not string or is empty.');
end

imname = deblank(imname);

% strip ,# volume numbers sometimes added by SPM
imname = regexprep(imname, ',\w*', '');

% strip .gz version for now and create .gz and non-.gz version. check for both.
imname = regexprep(imname, '.gz', '');

gzipexists = exist([imname '.gz'], 'file');
nozipexists = exist(imname, 'file');

if nozipexists
    % use imname
    
elseif gzipexists
    % use .gz - unzip and use previous value of imname
    % flag for re-gzipping later.
    waszipped = true;
    system(['gunzip ' imname '.gz']);
    
else
    error('Image "%s" does not exist or cannot be found on path (unzipped or .gz versions)!', imname);
end

full_imname = which(imname); % sometimes ok, sometimes empty if full path already included

if exist(full_imname, 'file') % if not empty and ok, else default to previous value of imname
    imname = full_imname;
end

% Add volume number if first vol only
if dofirstvolonly

    imname = [imname ',1'];
    
end

end
