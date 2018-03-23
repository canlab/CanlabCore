function [obj, obj_codes] = cat(obj, varargin)
% Merge two fmri_data objects, appending additional objects to the one
% entered first.  Resamples to the space of the first object if necessary.

% :Usage:
% ::
%
%     [obj, obj_codes] = cat(obj, [obj2], [obj3], etc.)
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2016 Tor Wager
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
%       An fmri_data object
%
%
% :Optional Inputs:
%   **obj2, obj3:**
%        A second/third/etc. fmri_data object to merge (concatenate) with the first
%   
% :Outputs:
%
%   **obj:**
%        The concatenated object
%
%   **obj_codes:**
%        Integer vector with one element for each image in the merged
%        object.  Integer values reflect which of the original objects the
%        image came from.
%
% :Examples:
% ::
%
%    If you have a series of fmri_data objects called DATA_OBJ in a cell array, 
%    concatenate them using: DATA_CAT = cat(DATA_OBJ{:});
%
%    Return image condition codes as well:
%    wh = [1 3];
%    [cat_obj, condition_codes] = cat(DATA_OBJ{wh});
%
% :References:
%   None
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    Created July 2016, Tor Wager
% ..


if isempty(varargin)
    
    obj_codes = ones(size(obj.dat, 2), 1);
    
    return
    
else
    
    obj_codes{1} = ones(size(obj.dat, 2), 1);
    
    for i = 1:length(varargin)
        
        obj_codes{i + 1} = (i + 1) .* ones(size(varargin{i}.dat, 2), 1);
        
        obj = merge_fcn(obj, varargin{i});
        
    end
    
    obj_codes = cat(1, obj_codes{:});
    
end

end % function






% SUBFUNCTIONS

function obj = merge_fcn(obj, obj2)

obj = replace_empty(obj,'voxels');   % changed to 'voxels' only. SG 2/23/18
obj2 = replace_empty(obj2,'voxels'); % changed to 'voxels' only. SG 2/23/18

obj2 = resample_space(obj2, obj);

% hcat these attributes
fnames = {'source_notes' 'Y_descrip' 'covariates_descrip' 'covariates' 'additional_info' 'dat'};

for i = 1:length(fnames)
    
    obj.(fnames{i}) = [obj.(fnames{i}) obj2.(fnames{i})];
    
end

% strvcat these attributes
fnames = {'image_names' 'fullpath'};

for i = 1:length(fnames)
    
    obj.(fnames{i}) = strvcat(obj.(fnames{i}), obj2.(fnames{i}));
    
end


% vcat these attributes
fnames = {'files_exist' 'X' 'Y' 'removed_images'};

for i = 1:length(fnames)
    
    obj.(fnames{i}) = [obj.(fnames{i}); obj2.(fnames{i})];
    
end

obj.history{end+1} = 'Concatenated another object with this one.';

end % function
