function [image_names, was_gzipped] = gunzip_image_names_if_gz(image_names)
% image_names = gunzip_image_names_if_gz(image_names)
%
% - Unzips images if they had .gz extension, and returns list of unzipped file names.
% - If unzipped file names are entered, returns names.
% - Returns string matrix or cell array, depending on input format (matches input)

% Tor Wager
% 5/5/17 Added deblank to extension
% 7/18/18 Added was_gzipped output for re-zipping

% Handle .gz by unzipping if needed

wascell = false;
was_gzipped = []; % do not pre-allocate because cell/char handled diff below

if iscell(image_names)
    wascell = true;
    image_names = char(image_names{:});
    
elseif ~ischar(image_names)
    
    return
    % not char array or cell of strings; nothing to do
    
end


image_names_cell = {};

for i = 1:size(image_names, 1)
    
    my_image = deblank(image_names(i, :));
    
    % Check if unzipped version already exists (name without .gz)
    unzipped_file = regexprep(my_image, '.gz', '');
    
    if exist(unzipped_file, 'file')
        %fprintf('Found zipped .gz file, but unzipped file \n%s \nalready exists. Using that file.\n', unzipped_file);
        
        image_names_cell{i, 1} = unzipped_file;
        
        % don't flag for re-zipping
        was_gzipped(i, 1) = false;
        
    else % We may have zipped version, or invalid file name
        
        [~, ~, ext] = fileparts(my_image);
        
        if strcmp(ext, '.gz')
            
            was_gzipped(i, 1) = true;
            
            %         image_names_cell(i, 1) = gunzip(my_image);
            
            % Unzip. Use system to remove .gz after unzipping.
            
            % Will wait for input, and not overwrite, if images exist
%             [status, result] = system(['gunzip ' my_image]);
            gunzip(my_image);

%             if status == 0 && isempty(result)
                % ok - success, and no error/text returned from system
                
                image_names_cell{i, 1} = unzipped_file;
                
%             else
%                 disp('Image unzipping error.')
%                 disp(result)
%             end
            
        else
            
            fprintf('Warning: Cannot find file:\n%s\n', my_image);
            
            was_gzipped(i, 1) = false;
            image_names_cell{i, 1} = my_image;
            
        end
    end
    
end % for each image

if wascell
    % do nothing
else
    image_names = char(image_names_cell{:});
end

end % function