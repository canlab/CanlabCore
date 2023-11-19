function [results_table_pos, results_table_neg, r, excluded_region_table, atlas_of_regions_covered, region_list, full_table] = table_of_atlas_regions_covered(obj, varargin)
% Make a table of which atlas parcels are covered by a set of regions or map identified in a study
%
% [results_table_pos, results_table_neg, r, excluded_region_table, atlas_of_regions_covered, region_list] = table_of_atlas_regions_covered(obj, [atlas_obj], ['percentile_threshold', 50])
%
% fMRI activation maps often include large suprathreshold areas that span
% multiple brain regions.  Traditional tables divide areas by contiguous
% regions ("blobs"), but this type of division is not useful if blobs cannot
% be easily summarized by a single anatomical label. A complementary approach
% is to identify which areas defined in a standard brain atlas ("parcels")
% are covered by the activation map. Parcels can be defined based on anatomy
% (e.g., cytoarchitecture) or function (e.g., resting-state fMRI).
% For example, a large blob may span the insula, claustrum, and putamen.
% This can be described in a table that lists all of these regions, along with
% the percentage of each parcel covered by the activation map.
%
% This function uses the object method subdivide_by_atlas() to create such
% a table, which is dividided into two tables with positive and negative average
% statistic values (results_table_pos, results_table_neg). These are Matlab
% table-class objects.  It also returns a region-class object for reference
% and rendering the subdivided blobs on brain slices or surfaces.
%
% You can enter any atlas-class object as a 2nd argument. If you don't,
% a default atlas object will be used.
%
% :Inputs:
%
%   **obj:**
%        A statistic_image,fmri_data, or image_vector object containing
%        only a single image to label. It should be thresholded for
%        meaningful results.
%
% :Optional Inputs:
%   **'coverage', coverage_percent_threshold:**
%        The keyword 'coverage' followed by a percentage threshold. 0 means
%        that any non-zero voxel in the input data that falls within an atlas
%        region will be sufficient to include the atlas region in the
%        table.  
%        50 means that 50% of the voxels in an atlas region must be
%        non-zero in the input image.
%        The default is 25%
%
%   **'atlas', [atlas_object]:**
%        An atlas-class object with regions
%
%   **noverbose, noprint:**
%        Suppress printing the tables
%
% :Outputs:
%
%   **results_table_pos, , r, excluded_region_table, table_legend_text:**
%        Table of results for regions with mean activation in map > 0
%
%   **results_table_neg:**
%        Table of results for regions with mean activation in map < 0
%
%   **r**
%        region-class object with covered regions
%
%   **excluded_region_table:**
%        Table of atlas regions below the coverage threshold
%
%   **atlas_of_regions_covered**
%        Atlas object with all regions covered
%
%   **region_list**
%        Text string with all covered regions and coverage percentages
%
%   **full_table**
%       Table of all regions, in order of indices in atlas       
%
% Examples:
% -----------------------------------------------------------------
% % Load a thresholded map:
% % (Note: You must have Neuroimaging_Pattern_Masks repo with subfolders on your path)
% % This is a predictive map for drug craving (Koban et al. 2022, Nat Neurosci)
%
% ncsthr = fmri_data(which('NCS_multithr_001k5_005_05_pruned.nii'), 'noverbose');
%
% [results_table_pos, results_table_neg, r, excluded_region_table, atlas_of_regions_covered, region_list, full_table] = ...
% table_of_atlas_regions_covered(ncsthr);
%
% [~, ~, r] = table_of_atlas_regions_covered(ncsthr, 'coverage', 50, 'noverbose');
%
% montage(r, 'regioncenters', 'colormap');
%
% Note: For a key to labels when using the Glasser 2016 Nature atlas, or the
% CANlab combined atlas (which uses Glasser), see GlasserTableS1.pdf in the
% original paper or CANlab github, or http://braininfo.rprc.washington.edu/
%
% see also region.table, for autolabeling of regions with an atlas and more
% detailed table output.

% ..
%     Author and copyright information:
%
%     Copyright (C) 2023 Tor Wager
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

% Tor Wager, Feb 2023
%            Nov 2023: Refactored to avoid using region method

if (size(obj.dat, 2) > 1)
    error('Use this function with fmri_data or image_vector objects containing only a single image. Use get_wh_image() to select one.')
end

% Old method:
% r = region(obj);
%
% [results_table_pos, results_table_neg, r, excluded_region_table, r_excluded, region_list] = table_of_atlas_regions_covered(r, varargin{:});

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% Defaults

percentile_threshold = 25;  % saves regions for which > this percentage of the atlas region is covered by the activation map
doprint = true;
atlas_obj = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'percentile_threshold', 'coverage', 'percent'}, percentile_threshold = varargin{i+1}; varargin{i+1} = [];

            case {'noverbose', 'nodisplay', 'noprint'}, doprint = false;

            case {'atlas_obj', 'atlas'}
                atlas_obj = varargin{i+1}; varargin{i+1} = [];

            case allowable_inputs
                % skip - handled above

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

validateattributes(percentile_threshold, {'double'}, {'scalar' 'nonnegative' 'nonnan'});

% -------------------------------------------------------------------------
% PREP ATLAS AND OBJECT
% -------------------------------------------------------------------------

% Load default atlas if needed
if isempty(atlas_obj)
    atlasname = 'canlab2023_combined_atlas_MNI152NLin6Asym_1mm.mat';
    atlaskeyword = 'canlab2023';

    if isempty(which(atlasname))

        %  atlasname = 'CANlab_combined_atlas_object_2018_2mm.mat';
        atlaskeyword = 'canlab2018_2mm';

    end

    atlas_obj = load_atlas(atlaskeyword);

end

validateattributes(atlas_obj, {'atlas'}, {});

% resample object to atlas
% this avoids dropping atlas regions and other complexities

obj = replace_empty(resample_space(obj, atlas_obj));
objdat = obj.dat;

if isa(obj, 'statistic_image')
    objdat = objdat .* obj.sig;    % apply current threshold
end

atlas_obj = replace_empty(atlas_obj);

n = length(atlas_obj.labels);

% -------------------------------------------------------------------------
% LOOP THROUGH REGIONS
% -------------------------------------------------------------------------

[atlas_region_coverage, n_atlas, mean_in_region, max_abs_in_region] = deal(zeros(n, 1));
regiondat = cell(1, n);

for i = 1:n

    indx = atlas_obj.dat == i;

    n_atlas(i, 1) = sum(indx);

    n_joint = sum(objdat & indx);

    % P(data | atlas) = n_joint / n_atlas

    atlas_region_coverage(i, 1) = 100 * n_joint / n_atlas(i);

    % mean and max in region
    regiondat{i} = obj.dat(objdat & indx);

    mean_in_region(i, 1) = nanmean(regiondat{i});

    max_abs_in_region(i, 1) = signedmax(regiondat{i});

    %     [mymax, whmax] = max(abs(regiondat{i}));
    %     if ~isempty(mymax)
    %         max_abs_in_region(i, 1) = sign(regiondat{i}(whmax)) .* mymax;
    %     end

end

atlas_region_coverage(n_atlas == 0) = 0;  % replace NaNs for convenience


% Build table
% ------------------------------------------------------

if ~iscolumn(atlas_obj.labels), atlas_obj.labels = atlas_obj.labels'; end

full_table = table(atlas_obj.labels, atlas_region_coverage, n_atlas, (1:n)',  mean_in_region, max_abs_in_region, 'VariableNames', {'Region' 'Coverage' 'Voxels_in_region' 'Atlas_index_number' 'mean_in_region' 'max_abs_in_region'});

% full_table = sortrows(full_table, {'Coverage' 'Voxels_in_region'}, 'descend');

output_table = full_table;
output_table(~(full_table.Coverage > percentile_threshold), :) = [];


% Get atlas regions
% ------------------------------------------------------
% note: this will reorder the numbers, in ascending order they were in the atlas
% so order of regions will not be the same. we need to reorder the table
% above, extract to match with the table, and re-order table and regions to
% sort by coverage again.

atlas_of_regions_covered = select_atlas_subset(atlas_obj, output_table.Atlas_index_number);

r = atlas2region(atlas_of_regions_covered);


output_table = sortrows(output_table, {'Atlas_index_number'}, 'ascend'); % to keep order-matched with region object, below
[output_table, indx] = sortrows(output_table, {'Coverage' 'Voxels_in_region'}, 'descend'); % indx maps from atlas number to coverage order

r = r(indx); % This should now be sorted into the same order as the table

% attach data values

for i = 1:length(r)

    r(i).Z = repmat(output_table.mean_in_region(i), 1, r(i).numVox);

end


% Separate tables and format output
% ------------------------------------------------------
results_table_pos = output_table(output_table.mean_in_region >= 0, :);

results_table_neg = output_table(output_table.mean_in_region < 0, :);

excluded_region_table = full_table(~(full_table.Coverage > percentile_threshold), :);

% r_excluded = [];  % not implemented


% Get region list
% -------------------------------------------------------------
region_list = cell(4, 1);
for i = 1:4, region_list{i} = ''; end
region_list{1} = 'Positive effects:';
region_list{3} = 'Negative effects:';

for i = 1:size(results_table_pos, 1)
    region_list{2} = [region_list{2} sprintf('%s (%2.0f%%), ', results_table_pos.Region{i}, results_table_pos.Coverage(i))];
end

for i = 1:size(results_table_neg, 1)
    region_list{4} = [region_list{4} sprintf('%s (%2.0f%%), ', results_table_neg.Region{i}, results_table_neg.Coverage(i))];
end

% Display table
% -------------------------------------------------------------
if doprint

    fprintf('Atlas name: %s\n', atlas_obj.atlas_name);
    fprintf('Coverage threshold: %3.1f%%\n\n', percentile_threshold);

    fprintf('\nPositive Effects\n\n')
    if ~isempty(results_table_pos)
        
        disp(results_table_pos)
    else
        disp('No regions to display');
    end

    fprintf('\nNegative Effects\n\n')
    if ~isempty(results_table_neg)
        
        disp(results_table_neg)
    else
        disp('No regions to display');
    end

end % doprint

end % main function




function val = signedmax(vals)

if isempty(vals), val = 0; return, end

vals = double(vals);
[maxabs, wh] = max(abs(vals));

val = sign(vals(wh)) .* maxabs;

end


