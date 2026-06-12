function obj = read_from_file(obj)
% read_from_file Read image data from disk into the .dat field of an image_vector object.
%
% Loads voxel data from the files listed in obj.fullpath into obj.dat.
% If obj.volInfo is empty, it is created from the first image. Try
% check_image_filenames(obj) first to make sure obj.fullpath points to
% valid files on disk.
%
% This is automatically called if you create a new image_vector object
% with names but do not directly enter data.
%
% :Usage:
% ::
%
%     obj = read_from_file(obj)
%
% :Inputs:
%
%   **obj:**
%        An image_vector / fmri_data object with .fullpath populated.
%        If .volInfo is empty, it will be created from the first image.
%
% :Outputs:
%
%   **obj:**
%        The input object with .dat populated as single-precision values
%        and .volInfo populated if it was empty.
%
% :Examples:
% ::
%
%     name = 'salientmap.nii';
%     img = image_vector('image_names', name);
%     % read_from_file is called automatically inside the constructor.
%
% :See also:
%   - check_image_filenames
%   - iimg_read_img
%   - iimg_get_data

% disp('Reading image data into object .dat field.');

% if we have mask info, use that. Otherwise, create
if isempty(obj.volInfo)
    % disp('Creating mask info from first image and storing in .volInfo.');
    obj.volInfo = iimg_read_img(obj.fullpath(1, :), 2);
else
    % disp('Using mask and space-defining info already stored in .volInfo field.');
end

% Now extract the actual data from the mask
switch spm('Ver')
    case {'SPM2', 'SPM99'}
        % legacy SPM
        obj.dat = iimg_get_data(obj.volInfo, obj.fullpath, 'single')';
    otherwise
        % SPM5+, including any future versions
        obj.dat = iimg_get_data(obj.volInfo, obj.fullpath, 'single', 'noexpand')';
end


end % function

