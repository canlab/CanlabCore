function write(obj, varargin)
%
% Write an image_vector object to an Analyze image.
% Option to write thresholded image, for statistic_image objects.
%
% obj.dat should contain data, with one COLUMN for each 3-D frame in the
% 4-D image to be written.
%
% Usage:
% write(obj)  -> writes to the image(s) specified in obj.fullpath
% write(obj, 'thresh') -> for statistic_image objects, writes thresholded
% write(obj, 'fname', '~/Documents/test.nii')  -> writes the image(s) to specific path
%
% For example: 
% If m is an image_vector object,
%
% m.X(m.X < .12) = 0; % apply an arbitrary but reasonable custom threshold
% orthviews(m);
% write the thresholded image to disk:
% anatmeanname = 'mean_gray_matter_mask.img';
% m.filename = anatmeanname;
% m.fullpath = fullfile(maskdir, anatmeanname);
% write(m)
%
% Option:
%   'mni'       resample image to standard MNI FOV (91x109x91)
%               uses mri_data.resample_space
%   'keepdt'    output image will be keep original data type (default = float32)
%   'fname'     writes out image to specific file name.  'fname' must be
%               followed by image name with path

%
% 2013/3/5: Luk[ea] added 'mni' option
% 2013/3/25: 
%           Luke[ea] added optional input to retain original datatype
% 2014/3/14: Luke added 'fname' option to specify filename

if any(strcmp(varargin, 'fname')) % fname option -- added by Luke
    obj.fullpath = varargin{find(strcmp(varargin, 'fname')) + 1}; %check if this works.
elseif isempty(obj.fullpath)
    error('Define fullpath field with name and path before writing an image_vector object to disk');
end

% Replace empty vox/images
obj = replace_empty(obj);

if any(strcmp(varargin, 'thresh'))
    if ~isa(obj, 'statistic_image')
        disp('Warning: Thresholding only works with statistic_image objects.');
    else
        disp('Writing thresholded statistic image.');
        obj.dat(~obj.sig) = 0;
    end
end

% mni option -- added Luk[ea]
if any(strcmp(varargin, 'mni'))
    evalc('mni = fmri_data(which(''brainmask.nii''));'); % evalc() used to silence output of fmri_data
    obj = resample_space(obj,mni);
end

% Keep original dt option -- added Luk[ea]
if any(strcmp(varargin,'keepdt'))
    iimg_reconstruct_vols(obj.dat, obj.volInfo, 'outname', obj.fullpath, 'keepdt');
else
    iimg_reconstruct_vols(obj.dat, obj.volInfo, 'outname', obj.fullpath);
end


fprintf('Writing: \n%s\n', obj.fullpath);



