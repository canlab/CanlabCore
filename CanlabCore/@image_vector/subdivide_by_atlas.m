% This function takes an activation map and subdivides it using an atlas
% passed in by the user. For example, if you have large clusters surviving
% thresholding, you can use this function to subdivide them along
% anatomical lines.
%
% :Usage:
% ::
%
%     [atlas_of_subdivisions, region_values] = subdivide_by_atlas(activation_map, atl, [optional inputs])
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2021 Tor Wager and Yoni Ashar
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
%   **activation_map:**
%        image to subdivide
%
%   **atlas:**
%        the subdivisions to apply to the activation_map
%
% :Optional Inputs:
%   **'table':**
%        NOT IMPLEMENTED YET
%        print a table of cluster coordinates by atlas regions
%
% :Outputs:
%
%   **atlas_of_subdivisions:**
%        an atlas; each atlas parcel is a cluster or part of a cluster from
%        the activation_map
%
%   **region_values:**
%        region object; each region corresponds to a parcel in
%        atlas_of_subdivisions and contains the activation values from 
%        the activation_map
%
% :Examples:
% ::
%
%    % give examples of code here
%    param1 = abc();
%    param2 = xyz();
%    [out1,out2] = func_call(param1, param2)
%
% :References:
%   CITATION(s) HERE
%
% :See also:
%   - list other functions related to this one, and alternatives*
%
function subdivided_atlas = subdivide_by_atlas(activation_map, atl)

    atl2 = apply_mask(activation_map, atl);
    r = atlas2region(atl2);
