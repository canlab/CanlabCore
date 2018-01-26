function [output_obj, voxel_matches] = select_voxels_by_value(parcel_obj, indx)
% Given an image_vector object and index values [integers], return an object
% with data values of 1 where voxel values match any index value, and 0
% otherwise. Esp. useful for selecting regions of interest from parcel/atlas images.
% 
% :Usage:
% ::
%
%     [output_obj, voxel_matches] = select_voxels_by_value(parcel_obj, vector of values)
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018 Tor Wager
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
%   **parcel_obj:**
%        an fmri_data object
%
%   **indx:**
%        a row vector of one or more values to search for, e.g., [5 12 17]
%
% :Optional Inputs:
%   None
%
% :Outputs:
%
%   **output_obj:**
%        Output object with 1/0 indicator vector as data
%
%   **voxel_matches:**
%        Matrix of which voxels match which index numbers, in the order you
%        entered the index numbers
%
% :Examples:
% ::
% Amygdala from Brainnetome atlas, Fan 2016
% parcellation_file = which('BN_Atlas_274_noCb_uint16.nii');
% parcel_obj = atlas(parcellation_file);
% output_obj = select_voxels_by_value(parcel_obj, [211:214]);
%
% [subregions, voxel_matches] = select_voxels_by_value(parcel_obj, [5 7 9]);
%
% :See also:
%   - fmri_data.apply_parcellation [method]
%

if ~isrow(indx), indx = indx'; end

objdat = double(parcel_obj.dat);

num_vox = size(objdat, 1);

indx_to_search = repmat(indx, num_vox, 1);

voxel_matches = objdat == indx_to_search;

% Find voxels
wh_vox = any(voxel_matches, 2);

output_obj = parcel_obj;
output_obj.dat = single(wh_vox);
output_obj = remove_empty(output_obj);


end % main function

