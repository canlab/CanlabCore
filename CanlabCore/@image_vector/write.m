function write(obj, varargin)
% Write an image_vector object to hard drive as an Analyze image (uses .fullpath field for image names)    
% Option to write thresholded image, for statistic_image objects.
%
% obj.dat should contain data, with one COLUMN for each 3-D frame in the
% 4-D image to be written.
%
% :Usage:
% ::
%
%    write(obj)  -> writes to the image(s) specified in obj.fullpath
%    write(obj, 'thresh') -> for statistic_image objects, writes thresholded
%    write(obj, 'fname', '~/Documents/test.nii')  -> writes the image(s) to specific path
%
% :Inputs:
%   **obj:**
%        An fmri_data, statistic_image, or other image_vector object.
%        Define a single character array with the path+file name to write
%        to in obj.fullpath.  It should end in .img or .nii to write a 3-D
%        or 4-D Analyze or Nifti file. If the file already exists, it will 
%        not be overwritten unless you also enter  the optional 'overwrite'
%        keyword.
%
% :Optional Inputs:
%
%   **mni:**
%        resample image to standard MNI mask dimensions (91x109x91, 2 mm vox)
%        uses mri_data.resample_space
%
%   **keepdt:**
%        output image will be keep original data type (default = float32)
%
%   **fname:**
%        writes out image to specific file name.  'fname' must be
%        followed by image name with path
%
%   **overwrite:**
%        Force overwrite of existing file. Use with caution!
%
% :Examples:
% ::
%
%    % If m is an image_vector object m.X(m.X < .12) = 0; % apply an 
%    % arbitrary but reasonable custom threshold
%    orthviews(m);
%
%    % write the thresholded image to disk:
%    anatmeanname = 'mean_gray_matter_mask.img';
%    m.filename = anatmeanname;
%    m.fullpath = fullfile(maskdir, anatmeanname);
%    write(m)
%

% ..
%    2013/3/5: Luk[ea] added 'mni' option
%
%    2013/3/25: Luke[ea] added optional input to retain original datatype
%
%    2014/3/14: Luke added 'fname' option to specify filename
%
%    2020/1: Tor changed to avoid evalc code. 
%
% ..

if any(strcmp(varargin, 'fname')) % fname option -- added by Luke
    
    obj.fullpath = varargin{find(strcmp(varargin, 'fname')) + 1}; 
    
elseif isempty(obj.fullpath)
    
    error('Define fullpath field with name and path before writing an image_vector object to disk');

end

% Check if file exists, and error if so, unless you use 'overwrite' option:

if exist(obj.fullpath, 'file') && ~any(strcmp(varargin, 'overwrite'))
    
    error('write() error: File already exists. Use ''overwrite'' option to force overwrite.');
    
end

% Check for illegal fullpath:

if size(obj.fullpath, 1) > 1
    
    error('Image name in obj.fullpath must be a single filename for 3-D or 4-D image.'); 
    
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
    
    % Edited 1/2020, evalc no longer needed
    %evalc('mni = fmri_data(which(''brainmask.nii''));'); % evalc() used to silence output of fmri_data
    mni = fmri_data(which('brainmask.nii'), 'noverbose');
    
    if isa(obj, 'statistic_image') % if obj is statistic_image, convert it to fmri_data first. 
        obj = fmri_data(obj);
        obj.mask = fmri_mask_image(obj);
    end
    
    obj = resample_space(obj, mni);
    
end

% Keep original dt option -- added Luk[ea]

if any(strcmp(varargin,'keepdt'))
    iimg_reconstruct_vols(obj.dat, obj.volInfo, 'outname', obj.fullpath, 'keepdt');
else
    iimg_reconstruct_vols(obj.dat, obj.volInfo, 'outname', obj.fullpath);
end


fprintf('Writing: \n%s\n', obj.fullpath);



