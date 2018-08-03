function [region_obj, region_table, table_legend_text] = autolabel_regions_using_atlas(region_obj, varargin)
% Load an atlas object (i.e., defined brain parcels with labels)  and use it to 
% label regions in a region object (e.g., 'blobs' or regions from an analysis)
%
% :Usage:
% ::
%
%     [region_obj, region_table] = autolabel_regions_using_atlas(region_obj, [atlas object with labels])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
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
%   **region_obj:**
%        A region object
%
% :Optional Inputs:
%   **atlas object with labels:**
%        An atlas-class object with labels to add. See @atlas, and
%        load_atlas
%
% :Outputs:
%
%   **region_obj:**
%        A copy of input region object with labels attached to the
%        .shorttitle and .title fields
%
%   **region_table:**
%        A table of summary statistics and more info on labels for each
%        region. region_table.Properties.Description and
%        region_table.Properties.VariableDescriptions have more info.
%
%   **table_legend_text:**
%        A descriptive legend with more info about the table.
%
% :Examples:
% -------------------------------------------------------------------------
% ::
% Example 1: 
% % Complete group analysis of a standard dataset
% % Do analysis and prep results region object:
%
%   img_obj = load_image_set('emotionreg');         % Load a dataset
%   t = ttest(img_obj, .005, 'unc');                % Do a group t-test
%   t = threshold(t, .005, 'unc', 'k', 10);         % Re-threshold with extent threshold of 10 contiguous voxels
%   r = region(t);                                  % Turn t-map into a region object with one element per contig region
%   montage(r);                                     % show slices; see region.montage
%
%   Label regions and print a table:
%   [r, region_table, table_legend_text] =
%   autolabel_regions_using_atlas(r);               % Label regions. Can be skipped because 'table' below attempts to do this automatically
%   table(r);                                       % Print a table of results using new region names
%
% :References:
%   See table_legend_text and atlas object loaded for references for published atlases used.
%
% :See also:
%   - atlas_similarity.m (main function used here, for atlas objects)
%   - table method for region object (uses this function)
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..
% Test code: Validate by visual inspection of each region and label
% 
% for i = 1:length(region_obj)
%     
%     orthviews(region_obj(i), {[1 0 0]});
%     
%     myindx = region_table.modal_atlas_index(i);
%     if myindx
%         test = select_atlas_subset(atlas_obj, myindx);
%         orthviews(test, 'color', {[0 0 1]}, 'add', 'largest_cluster');
%         disp(region_table.modal_label{i})
%     else
%         disp('No match')
%     end
%     
%     a = input('press');
% end


% Load atlas
% ------------------------------

if nargin < 2
    
 % 2 mm loads MUCH faster than 1 mm
atlas_obj = load_atlas('canlab2018_2mm', 'noverbose'); 

else
    
    atlas_obj = varargin{1};
    
end

% Prep atlases and match spaces
% ------------------------------

r = region2atlas(region_obj);

k_orig = length(region_obj); % original regions

% Resample atlas space so they match
% Will compress indices and remove regions that do not match
% atlas_obj = resample_space(atlas_obj, r); % resample atlas to region space - really slow if atlas_obj is higher-res

r = resample_space(r, atlas_obj); 

k = num_regions(r); % after resampling

% Check for loss of regions -- this could create errors in labeling:

if k ~= k_orig
    disp('Regions lost due to resampling. Check and develop code.');
    keyboard
end

% Main work done here, on atlases
% ------------------------------
[region_table, table_legend_text, coverage25_labels] = atlas_similarity(r, atlas_obj);

% Add region labels
% ------------------------------

for i = 1:length(region_obj)
    
    region_obj(i).shorttitle = region_table.modal_label{i};
    
    region_obj(i).title = char(coverage25_labels{i});
    
end % function
