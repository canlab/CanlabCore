function obj_selected = select_regions_near_crosshairs(atlas_obj, varargin)
% atlas object method: Select regions within x mm of the crosshairs of an spm_orthviews display
% :Usage:
% ::
%
%     obj_selected = select_regions_near_crosshairs(obj, [distance_threshold])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2019 Tor Wager
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
%        an atlas object
%
% :Optional Inputs:
%   **thresh:**
%        'thresh' followed by a distance threshold in mm. Default is 15 mm
%
%   **coords:**
%        'coords' followed by an [x y z] coordinate triplet in mm
%
% :Outputs:
%
%   **obj_selected:**
%        an atlas object containing the nearby regions
%
% :Examples:
% ::
%
%    % Load an atlas and find coordinates within 10 mm
%    % -------------------------------------------------------------
%    atlas_obj = load_atlas('canlab2018_2mm');
%    orthviews(atlas_obj)
%    obj_within_10mm = select_regions_near_crosshairs(atlas_obj, 'thresh', 10)
%    obj_within_10mm.labels
%
%    % Load an atlas and find coordinates within 15 mm of [-40 0 20]
%    % -------------------------------------------------------------
%    atlas_obj = load_atlas('canlab2018_2mm');
%    orthviews(atlas_obj)
%    obj_within_15mm = select_regions_near_crosshairs(atlas_obj, 'coords', [-40 0 20])
%    obj_within_15mm.labels
%    
%
% :References:
%   None.
%
% :See also:
%   - cluster_graphic_select, for a non-object-oriented selector
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% ..
%    DEFAULTS AND INPUTS
% ..

% Some useful things
thresh = 15;
coords = [];

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'thresh', thresh = varargin{i+1}; varargin{i+1} = [];
            case 'coords', coords = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(coords)
    coords = spm_orthviews('pos')';
end

obj = atlas2region(atlas_obj); % convert to region
centers = cat(1,obj.mm_center);

% find  cluster, based on center
d = distance_euclid(coords,centers); 

wh = find(d <= thresh);

obj_selected = select_atlas_subset(atlas_obj, wh);

end

