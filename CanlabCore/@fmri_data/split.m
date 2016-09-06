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

en = cumsum(imgs_per_obj);
st = [0 en(1:end-1)] + 1;

k = length(imgs_per_obj);
obj_cell = cell(1, k);

for i = 1:k
    
    obj_cell{i} = obj;
    
    wh = st(i):en(i);
    
    obj_cell{i}.dat = obj_cell{i}.dat(:, wh);
    
    obj_cell{i}.images_per_session = en(i) - st(i) + 1;
    
    if ~isempty(obj_cell{i}.image_names) && size(obj_cell{i}.image_names, 1) > en(i)
        
        obj_cell{i}.image_names = obj_cell{i}.image_names(wh, :);
        
    end
    
    if ~isempty(obj_cell{i}.fullpath) && size(obj_cell{i}.fullpath, 1) > en(i)
        
        obj_cell{i}.fullpath = obj_cell{i}.fullpath(wh, :);
        
    end
    
    if ~isempty(obj_cell{i}.files_exist) && size(obj_cell{i}.files_exist, 1) > en(i)
        
        obj_cell{i}.files_exist = obj_cell{i}.files_exist(wh, :);
        
    end
    
    obj_cell{i} = check_image_filenames(obj_cell{i}, 'noverbose');
    
end

end % function
