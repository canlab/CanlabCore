function [r, local_pattern_response] = extract_data(r, data_obj)
% Extract data and apply local patterns from image_vector object (data_obj) for voxels specified by a region object (r). Returns extracted data and averages.
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
%    - data_obj does not have to be in the same space, uses mm coordinates
%    - if 2nd output is requested and r.vals is nonempty, will compute
%    local pattern responses for each region, using r(i).vals as weights.
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
%   **local_pattern_response:**
%         local linear combinations of voxels in each region, computed
%         using r(i).vals as weights.
%
% ..
%     Programmers' notes:
%     8/3/2015 Tor Wager: Fixed bug in averaging when only 1 voxel in region
%     5/15/2017 Tor Wager : Added replace_empty to avoid voxel list
%     mismatch when empty vox were removed
% ..

if isempty(r), return, end

xyzlist = data_obj.volInfo.xyzlist;

data_obj = replace_empty(data_obj); % tor added - May 2017

k = size(data_obj.dat, 2);  % images to extract from

for i = 1:length(r)
    
    % get region voxel coords, but in space of data_obj
    regionxyz = mm2voxel(r(i).XYZmm, data_obj.volInfo.mat, 1); % allows repeats
    
    v = size(regionxyz, 1); % num voxels
    
    % find wh_vox, list of voxels in data object that match region coords
    [whregionxyz, wh_vox] = match_coordinates(regionxyz, xyzlist); 

    % build data, using NaN where we cannot find a match
    dat = NaN .* zeros(k, v);
    dat(:, whregionxyz) = data_obj.dat(wh_vox, :)';
    
    weight_vals = double(r(i).val);      % linear combo / pattern expression
    
    r(i).all_data = dat;
    r(i).dat = nanmean(dat, 2);
    r(i).source_images = data_obj.fullpath;
    
    % linear combo / pattern expression
    % -----------------------------------------
    if nargout > 1 && ~isempty(weight_vals)
    % then apply r(i).vals ...
    
        go_ok = true;
        if length(weight_vals) ~= v
            disp('Warning: weights in r(i).dat not empty, but wrong length.');
            go_ok = false;
        end
        
        whnan = isnan(dat(:));
        if any(whnan)
            sprintf('Warning: missing voxels in region %3.0f. Applying pattern to existing data\n', i);
            dat(whnan) = 0;
        end
        
        whnan = isnan(weight_vals);
        if any(whnan)
            sprintf('Warning: missing voxels in weight vector for region %3.0f. Applying existing weights\n', i);
            weight_vals(whnan) = 0;
        end
        
        if go_ok
            local_pattern_response{i} = double(dat) * weight_vals;
        end
        
    end
end

end % function

function  [whregionxyz, wh_vox] = match_coordinates(regionxyz, xyzlist)

    % wh_vox is index into data_obj
    % whregionxyz is index into region voxels
    
    [~, whregionxyz, wh_vox] = intersect(regionxyz, xyzlist, 'rows');

end
