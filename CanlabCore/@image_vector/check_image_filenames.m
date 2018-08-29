function obj = check_image_filenames(obj, varargin)
% Check whether images listed in obj.fullpath actually exist
%
% :Usage:
% ::
%
%     obj = check_image_filenames(obj, ['noverbose'])
%
% :Behavior:
%    - If there are no file names, do nothing.
%    - If file names are entered and full path is not, attempt to find full
%      path.
%    - If full path info is entered, check to see if files exist.
%      Return output in obj.files_exist, and print a warning if only some exist.
%
% Image names should be stored in .fullpath
% abbreviated image names may be stored in image_names.
%
% :Note:
% 
% fullpath should have full path to each volume in a string matrix, with
% trailing [,volume#] for 4-D images as per SPM style expanded list.
%
% image_names should have image name only for each volume
%
% ..
%    May still be debugging issues with 3-D vs. 4-D files
% ..
%
% Tor Wager: 7/2018 Added support for .gz file checking, avoiding spm_vol
% warning

verbose = isempty(strmatch('noverbose', varargin(cellfun(@ischar, varargin)))); % if 'noverbose' is entered, suppress output

% Make char array if needed
if ~isempty(obj.image_names) && iscell(obj.image_names)
    
    obj.image_names = char(obj.image_names{:});
    
end

if isempty(obj.fullpath) && ~isempty(obj.image_names)
    % get full path
    if verbose, disp('Image names entered, but fullpath attribute is empty. Getting path info.'); end
    fullnames = [];
    imgnames = [];
    
    for i = 1:size(obj.image_names, 1)
        
        if any(obj.image_names(i, :) == filesep)
            
            % may already have full path
            fullnames = strvcat(fullnames, deblank(obj.image_names(i, :)));
            [dd, ff, ee] = fileparts(deblank(obj.image_names(i, :)));
            imgnames = strvcat(imgnames, [ff ee]);
            
        else
            thisname = deblank(obj.image_names(i, :));
            fullnames = strvcat(fullnames, which(thisname));
            %             if exist(which(thisname), 'file')
            %             elseif exist(fullfile(pwd, thisname), 'file')
            %                 fullnames = strvcat(fullnames, fullfile(pwd, thisname));
            %             end
            imgnames = strvcat(imgnames, thisname);
        end
    end
    
    obj.fullpath = fullnames;
    obj.image_names = imgnames;
    
    %     disp('Warning: image names with full paths should be stored in obj.fullpath.')
    %     disp('Your fullpath field is empty, but you have some names in image_names.');
    %     disp('These names will not be checked, and this program assumes the images');
    %     disp('corresponding to your dataset do not exist.');
end

if isempty(obj.fullpath)
    obj.files_exist = false(max([1 size(obj.image_names, 1) size(obj.fullpath, 1)]), 1);
end

% Check if same size as number of volumes in data
if ~isempty(obj.fullpath) && size(obj.fullpath, 1) ~= size(obj.dat, 2)
    
    if verbose
        disp('.fullpath should have image name for each image column in .dat');
        disp('Attempting to expand image filenames in case image list is unexpanded 4-D images');
    end

    if any(strfind(obj.fullpath(1,:), '.gz'))
        % skip expand_4d_filenames, cause of spm_vol (returns a warning, but works...)
        % instead, list filenames manually.
        obj.fullpath = repmat(obj.fullpath, size(obj.dat, 2), 1);
    else
        
        obj.fullpath = expand_4d_filenames(obj.fullpath);
        
    end
    
    
end

for i = 1:size(obj.fullpath, 1)
    
    img = deblank(obj.fullpath(i, :));
    
    % take off trailing commas (SPM notation for multiple vols)
    wh = find(img == ',');
    if ~isempty(wh)
        wh = wh(end);
        img(wh:end) = [];
    end
    
    obj.files_exist(i, 1) = exist(img, 'file') > 0;
    
end

% If some, but not all, files exist, print a warning:
if any(obj.files_exist)
    if sum(obj.files_exist) < size(obj.fullpath, 1)
        if verbose
            disp('Warning: *Some* files are missing. Printing index numbers of missing files:');
            disp(find(~obj.files_exist));
        end
    end
end

obj.history{end+1} = sprintf('Checked image files exist. All exist = %3.0f', ...
    all(obj.files_exist));

end % function


