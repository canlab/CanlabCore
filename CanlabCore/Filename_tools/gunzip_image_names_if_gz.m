function image_names = gunzip_image_names_if_gz(image_names)
% image_names = gunzip_image_names_if_gz(image_names)
%
% - Unzips images if they had .gz extension, and returns list of unzipped
% file names.
% - Returns string matrix or cell array, depending on input format (matches
% input)

% Handle .gz by unzipping if needed

wascell = false;

if iscell(image_names)
    wascell = true;
    image_names = char(image_names{:});
    
elseif ~ischar(image_names)
    return
    % not char array or cell of strings; nothing to do
end
    

image_names_cell = {};

for i = 1:size(image_names, 1)
    
    [~, ~, ext] = fileparts(image_names(i, :));
    
    if strcmp(ext, '.gz')
        image_names_cell(i, 1) = gunzip(image_names(i, :));
    else
        image_names_cell{i, 1} = image_names(i, :);
    end
end

if wascell
    % do nothing
else
    image_names = char(image_names_cell{:});
end

end % function