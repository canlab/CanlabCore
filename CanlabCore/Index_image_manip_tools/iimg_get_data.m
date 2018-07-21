% Read data within a mask (described by maskInfo) for 3-D or 4-D image data
% See examples below.
%
% :Usage:
% ::
%
%     [dat, maskInfo] = iimg_get_data(mask, imageNames, varargin)
%
% :Inputs:
%
%   3 input options for mask:
%      1) a filename string, e.g., 'mask.img'
%      2) volInfo structure (from iimg_read_img)
%      3) xyz list, n rows x 3 voxel coordinates
%
% :Optional inputs:
%
%   **'nocheck':**
%        Skip checking of image dimensions for homogeneity (saves time)
%
%   **'verbose':**
%        Verbose output
%
%   **'single':**
%        This option may be helpful for large datasets.  It
%        pre-allocates a single-precision matrix for the data (which
%        saves space), and then loads the data one image at a time.
%
%   **'noexpand':**
%        Skip expanding the list of volumes for 4-D filenames.  The
%        expansion is the default, and is compatible with spm2 and spm5/8,
%        but is not needed for spm5 and above.
% 
% :Examples:
% ::
%
%    [dat, volInfo] = iimg_get_data('graymask.img', imgs);
%
%    xyz = [20 20 20; 25 25 25; 30 30 30; 5 5 5];
%    [dat, volInfo] = iimg_get_data(xyz, imgs);
%
% Example of image reading
%
%   1) Get volume info from first volume of a 3-D or 4-D image
%      volInfo structure has necessary info for converting to/from
%      "vectorized" format
%      dat returns 4-D data from entire image (all volumes)
%      ::
%
%          img_name = 'test_run1_pca.img';
%          [maskInfo, dat] = iimg_read_img(img_name, 2);
%
%   2) Get data in "vectorized" image format for each volume in the
%      image. Works for a list of images too. data is in-mask voxels x volumes
%      this can be useful when you want to return whole-brain data for many
%      images, but in a search volume only (i.e., no extra-brain voxels)
%      ::
%
%          data = iimg_get_data(maskInfo, img_name);
%          data = data';  % make sure columns are volumes
%
%   3) Write out a 4-D image with the same data, called test_run1_pca2.img
%      ::
%
%          voldat3D = iimg_reconstruct_vols(data, maskInfo, 'outname',
%                                                  'test_run1_pca2.img');

function [dat, maskInfo] = iimg_get_data(mask, imageNames, varargin)
docheck = 1;
verbose = 0;
dosingle = 0;
doexpand = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'nocheck', docheck = 0;
            case 'verbose', verbose = 1;
            case 'noverbose'             % nothing more to do
            case 'single', dosingle = 1; % varargin{i+1};
            case 'noexpand', doexpand = 0;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% get mask info structure and image info
% ------------------------------------------------------------
if verbose, fprintf('loading mask. '); tic, end

if isstruct(mask)
    maskInfo = mask;
    if ~isfield(maskInfo, 'xyzlist')
        error('Need xyzlist in maskInfo.  Try using extended output flag with iimg_read_img.');
    end
elseif ischar(mask)
    maskInfo = iimg_read_img(mask, 1);
else
    % xyz list
    maskInfo = struct('xyzlist', mask);
    docheck = 0;
end

wasempty = isempty(imageNames);

% Expand 4-D filenames for both SPM2 and SPM5 compatibility
if doexpand
    if verbose, fprintf('expanding image name list. '); end
    imageNames = expand_4d_filenames(char(imageNames)); % for SPM2 compatibility
end

if isempty(imageNames)
    if wasempty, error('Image name list is empty');
    else error('Image names do not seem to exist on path');
    end
end

if verbose, fprintf('mapping volumes. '); end
imgInfo = spm_vol(imageNames);

% Check the dimensions
% ------------------------------------------------------------
% 1/14/17: Check ALL images, in case we have stacked images in different
% spaces into the same list
if docheck
    if verbose, fprintf('\nchecking that dimensions and voxel sizes of volumes are the same. '); end
    % anybad = iimg_check_volinfo(maskInfo, imgInfo(1));
    anybad = iimg_check_volinfo(maskInfo, imgInfo);
    
    if anybad
        disp('Reslice images so dimensions and vox sizes match.');
        disp('Try using the function scn_map_image for an easy way to do this.');
        error('Exiting.');
    end
end

if dosingle
    nvols = length(imgInfo);
    nvox = size(maskInfo.xyzlist, 1);
    if verbose, fprintf('\nPre-allocating data array. Needed: %3.0f bytes\n', 4*nvols*nvox); end
    dat = zeros(nvols, nvox, 'single');

    % load one-at-a-time into single matrix, to save memory
    if verbose, fprintf('Loading image number: %5.0f', 0); end

    for i = 1:nvols
        if verbose, fprintf('\b\b\b\b\b%5.0f', i); end
        dat(i, :) = spm_get_data(imgInfo(i), maskInfo.xyzlist');
    end

    if verbose, fprintf('\n'); end
else
    if verbose, fprintf('\nLoading data into double-precision array'); end
    dat = spm_get_data(imgInfo, maskInfo.xyzlist');

end

if verbose, toc, end
% Notes:
% This is faster for lists with many voxels (like 200,000)
% but is slower for 44,000 voxels
%tic, dat2 = spm_read_vols(spm_vol(imageNames));, toc
end
