function obj_cell = split(obj, varargin)
% Split fmri_data object into separate objects with different images in each
% default: Use the obj.images_per_session field
%
% :Usage:
% ::
%
%     obj_cell = split(obj, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C)2016 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **obj:**
%        An fmri_display object with multiple images, to split
%
% :Optional Inputs:
%   **'imgs_per_obj'**
%        Followed by a vector of number of images in each object
%        e.g., for 100 images, [25 25 25 25] would split into 4 objects of
%        25 images each
%        Default: Use obj.images_per_session
%
% :Outputs:
%
%   **objcell:**
%        A cell array with each output fmri_data object in a cell.
%
%
% :Examples:
% ::
%
%
% :References:
%   None.
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..
% Tor Wager, Aug 2016
% Added  >= en(i) to some lines - bug fix. Seemed to have little effect on most cases. Tor, June 2017

% ..
%    DEFAULTS AND INPUTS
% ..

imgs_per_obj = obj.images_per_session;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'imgs_per_obj', imgs_per_obj = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(imgs_per_obj)
    obj_cell = {obj};
    return
end

obj = replace_empty(obj);  % to replace removed images

if ~isrow(imgs_per_obj), imgs_per_obj = imgs_per_obj'; end

en = cumsum(imgs_per_obj);
st = [0 en(1:end-1)] + 1;

k = length(imgs_per_obj);
obj_cell = cell(1, k);

% warn user about fields that aren't handled but appear like they should be
% split. Some of these (e.g. additional_info) are handled by
% fmri_data.cat(), so it may surprise the user that split doesn't invert
% the process.
%
% criteria for warning: if length of field is equal to number of splits
fnames = {'image_names', 'fullpath', 'files_exist', 'X', 'Y'};
supportedFields = {fnames{:}, 'dat', 'images_per_session', 'removed_images'};
allFields = fieldnames(obj);
for i = 1:length(allFields)
    if sum(strcmp(allFields{i},supportedFields)) == 0
        if length(obj.(allFields{i})) == k
            warning(['Object''s ''', allFields{i}, ''' field length matches number of splits, but ''', allFields{i}, ''' field will not be split'])
        end
    end
end
    
for i = 1:k
    
    obj_cell{i} = obj;
    
    wh = st(i):en(i);
    
    obj_cell{i}.dat = obj_cell{i}.dat(:, wh);
    
    obj_cell{i}.images_per_session = en(i) - st(i) + 1;
    
    
    for j=1:length(fnames)

        if ~isempty(obj_cell{i}.(fnames{j})) && size(obj_cell{i}.(fnames{j}), 1) >= en(i)

            obj_cell{i}.(fnames{j}) = obj_cell{i}.(fnames{j})(wh, :);

        end
    end
    
    obj_cell{i} = check_image_filenames(obj_cell{i}, 'noverbose');
    

    if ~isempty(obj_cell{i}.removed_images) && size(obj_cell{i}.removed_images, 1) >= en(i)
        
        obj_cell{i}.removed_images = obj_cell{i}.removed_images(wh, :);
        
    else

        obj_cell{i}.removed_images = 0; % add in case

    end    

end

end % function
