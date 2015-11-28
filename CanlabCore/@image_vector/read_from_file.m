function obj = read_from_file(obj)
% Reads data from image filenames into obj.dat
%
% Try obj = check_image_filenames(obj) first.
%
% This is automatically called if you create a new image_vector object with
% names but do not directly enter data. e.g., the commands below will load data:
%   - name = 'salientmap.nii';
%   - img = image_vector('image_names', name);
% 

disp('Reading image data into object .dat field.');

% if we have mask info, use that. Otherwise, create
if isempty(obj.volInfo)
    disp('Creating mask info from first image and storing in .volInfo.');
    obj.volInfo = iimg_read_img(obj.fullpath(1, :), 2);
else
    disp('Using mask and space-defining info already stored in .volInfo field.');
end

% Now extract the actual data from the mask
switch spm('Ver')
    
    case {'SPM8', 'SPM5','SPM12'}
        obj.dat = iimg_get_data(obj.volInfo, obj.fullpath, 'single', 'noexpand')';
        
    case {'SPM2', 'SPM99'}
        % legacy, for old SPM
        obj.dat = iimg_get_data(obj.volInfo, obj.fullpath, 'single')';
        
    otherwise
        error('Unknown version of SPM! Update code, check path, etc.');
end


end % function

