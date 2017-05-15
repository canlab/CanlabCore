function r = extract_data(r, data_obj)
% Extract data from image_vector object (data_obj) for voxels specified
% by a region object (r). Returns extracted data and averages.
%
%
% :Usage:
% ::
%
%    region_obj = extract_data(region_obj, data_obj)
%
% Type methods(region) for a list of special commands for region object
% Type help object_name.method_name for help on specific methods.
%
% :Features:
%    data_obj does not have to be in the same space, uses mm coordinates
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2010 Tor Wager
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
%   **r:**
%         a region object
%
%   **data_obj:**
%         an image_vector or fmri_data object to extract data from
%         does not have to be in the same space, uses mm coordinates
%         NOTE: resampling the image first will give slightly different
%         results, due to interpolation
%
% :Outputs:
%
%   **r:**
%         a region object, with data attached
%
% ..
%     Programmers' notes:
%     8/3/2015 Tor Wager: Fixed bug in averaging when only 1 voxel in region
%     5/15/2017 Tor Wager : Added replace_empty to avoid voxel list
%     mismatch when empty vox were removed
% ..

xyzlist = data_obj.volInfo.xyzlist;

data_obj = replace_empty(data_obj); % tor added - May 2017

k = size(data_obj.dat, 2);  % images to extract from

for i = 1:length(r)
    
    % get region voxel coords, but in space of data_obj
    regionxyz = mm2voxel(r(i).XYZmm, data_obj.volInfo.mat, 1); % allows repeats
    
    v = size(regionxyz, 1); % num voxels
    
    [whregionxyz, wh_vox] = match_coordinates(regionxyz, xyzlist); 

    % build data, using NaN where we cannot find a match
    dat = NaN .* zeros(k, v);
    dat(:, whregionxyz) = data_obj.dat(wh_vox, :)';
    
    r(i).all_data = dat;
    r(i).dat = nanmean(dat, 2);
    r(i).source_images = data_obj.fullpath;
end

end % function

function  [whregionxyz, wh_vox] = match_coordinates(regionxyz, xyzlist)

    % wh_vox is index into data_obj
    % whregionxyz is index into region voxels
    
    [~, whregionxyz, wh_vox] = intersect(regionxyz, xyzlist, 'rows');

end
