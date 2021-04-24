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
%   **obj:**
%        image to subdivide 
%        e.g., image_vector, fmri_data or statistic_image object containing activation map
%        single-image objects only. See get_wh_image( ) to select one from a set.
%
%   **atlas:**
%        atlas object with the subdivisions to apply to the activation_map
%        see load_atlas( ) for examples/pre-defined atlases.
%
% :Optional Inputs:
%   **'table':**
%        NOT IMPLEMENTED YET
%        print a table of cluster coordinates by atlas regions
%
% :Outputs:
%
%   **atlas_of_subdivisions:**
%        an atlas object
%        subdivided_atlas contains voxel in the original map AND the atlas, partitioned by atlas regions.
%        each atlas parcel is a cluster or part of a cluster from the activation_map
%
%   **region_values:**
%        region object; each region corresponds to a parcel in
%        atlas_of_subdivisions and contains the activation values from 
%        the activation_map
%
% :Examples:
% ::
%
% % Example 1: Subdivide a thresholded t-statistic image
% % --------------------------------------------------------------------
% Load sample images, creating and fmri_data class object with 30 images     
% 	imgs = load_image_set('emotionreg');
% 
% % Perform a t-test on each voxel, returning a statistic_image object
% % containing t-stats and p-values:
%	t = ttest(imgs);
% 
% % Threshold the t-statistic_image object at p < 0.005
%	t = threshold(t, .005, 'unc');
% 
% % Load an atlas object whose boundaries we want to apply
%
%   atl = load_atlas('painpathways');
% 
% % Subdivide the map into regions defined by the atlas.
% % subdivided_atlas contains voxel in the original map AND the atlas, partitioned by atlas regions.
%
%   [subdivided_atlas, r] = t.subdivide_by_atlas(atl);
%
% % Plot the results
% 
%   figure; montage(subdivided_atlas);
%   montage(r, 'regioncenters', 'colormap');
%
% % Make a table of results, with one row per atlas-defined region:
%
%   [results_table, results_table_pos, results_table_neg] = table_simple(r);
% 
% % Make a more elaborate table labeled by a potentially different atlas:
% % Note: the standard table may use a DIFFERENT atlas, chosen in the table( )
% method. So the "atlas regions covered" here may refer to a different atlas.
%
%   [poscl, negcl, results_table] = table(r, 'nolegend');
%
% % Example 2: Subdivide a mask image
% % --------------------------------------------------------------------
% % Load an activation map or mask (e.g., thresholded results image)
%
% img = '/Users/torwager/Documents/GitHub/OLP4CBP/data/bladder_pain_param_mod/evoked_pain_localizer_masked.nii';
% obj = fmri_data(img);
% atl = load_atlas('painpathways');
% [subdivided_atlas, r] = subdivide_by_atlas(obj, atl);
% [results_table, results_table_pos, results_table_neg] = table_simple(r);
%
% :References:
%   N/A
%
% :See also:
%   - region.subdivide_by_atlas
%
function [subdivided_atlas, r] = subdivide_by_atlas(obj, atl)

% Check inputs
% ---------------------------------------------

% Error if > 1 image
if size(obj.dat, 2) > 1
    error('Use image_vector.subdivide_by_atlas with single-image objects only. See get_wh_image( ) to select.');
end

% Main work done here
% ---------------------------------------------

subdivided_atlas = apply_mask(atl, obj);
r = atlas2region(subdivided_atlas);

% Extract the original image values from the activation map and save them
% in region object, for table-making, etc.
% ---------------------------------------------

% Some voxels are in obj but not mask; exclude these
obj2 = resample_space(obj, subdivided_atlas);
obj2 = apply_mask(obj2, subdivided_atlas);

% Make sure we have the same voxel list, squeezing empties from both
subdivided_atlas = remove_empty(subdivided_atlas); 
obj2 = remove_empty(obj2);

for i = 1:length(r)

    wh = subdivided_atlas.dat == i;
    
    vox_vals = obj2.dat(wh, 1);
    
    % Replace .val and .Z both, because they are/may be used in
    % table-making, and they're (unfortunately) redundant. Z is legacy.
    
    r(i).val = vox_vals;
    r(i).Z = vox_vals';
    
end

end % main function

