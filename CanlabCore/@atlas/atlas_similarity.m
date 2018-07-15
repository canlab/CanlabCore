function [region_table, coverage50_labels, coverage50_index, x_counts, x_dice, x_atlas_coverage] = atlas_similarity(atlas_to_parse, ref_atlas_obj)
% Take regions in an atlas object (atlas_to_parse) and annotate them with labels and quantitative
% coverage stats from another atlas (ref_atlas_obj)
%
% [region_table, coverage50_labels, coverage50_index, x_counts, x_dice, x_atlas_coverage] = atlas_similarity(atlas_to_parse, ref_atlas_obj)
%
% Dice: Compute similarities (cross-counts and Dice coeffs) between parcels in two atlases
% matrix of [atlas1] x [atlas2]
%
% Mode: Labels each region in a target atlas (atlas 1) with best-matching (modal) label
% in a reference atlas (atlas 2). 
% - modal_atlas_index is the reference atlas index number for each target region
% - modal_atlas_label is the reference atlas label for each target region, or "No
% region identified" if the target region does not match any reference
% atlas region.
%
% Coverage: Return matrix of coverage of reference atlas regions for each
% target region. For each target (row), values are the proportion of
% reference region (column) covered. 
% - modal_atlas_coverage is the proportion of reference atlas voxels
% covered by the best-matching region
%
% Labels:

% Percent of voxels in each atlas region covered by the blob
% Intersection / size of atlas region

% note: vox counts can be zero if regions are (1) outside mask, or (2)
% empty after reslicing to the atlas space

% note: omits 0
[~, voxcount_atlas, imtx_atlas] = get_region_volumes(ref_atlas_obj);

[vol_regions, voxcount_regions, imtx_regions] = get_region_volumes(atlas_to_parse);

% Cross-counts
% ------------------------------------------------------

x_counts = imtx_regions' * imtx_atlas;  % counts of regions x atlas parcels

% Dice coefficients
% ------------------------------------------------------
nr = length(voxcount_regions);
na = length(voxcount_atlas);

sum_for_region = repmat(voxcount_regions', 1, na);
sum_for_atlas = repmat(voxcount_atlas, nr, 1);

x_dice = x_counts ./ (sum_for_region + sum_for_atlas);

% Coverage of each atlas region
% ------------------------------------------------------

x_atlas_coverage = x_counts ./ sum_for_atlas;


% Modal atlas region for each input region
% ------------------------------------------------------

[max_count_by_region, modal_atlas_index] = max(x_counts, [], 2);

wh_no_atlas_region = max_count_by_region == 0;

modal_atlas_index(wh_no_atlas_region) = 0;  % because max returns 1 where empty counts 

% Modal labels
wh_labels = modal_atlas_index > 0;

modal_label = repmat({'No region identified.'}, nr, 1);
modal_label(wh_labels) = ref_atlas_obj.labels(modal_atlas_index(wh_labels));

% modal_atlas_coverage: Percent of best reference region covered
%
for i = 1:nr
   
    if wh_no_atlas_region(i) 
        % Missing - no atlas regions match
        modal_atlas_coverage(i, 1) = 0;
    else
        mycoverage = x_atlas_coverage(i, :);
        modal_atlas_coverage(i, 1) = mycoverage(modal_atlas_index(i));
    end
    
end

% Percent in modal atlas region


% Reference regions with at least 50% coverage
% ---------------------------------------------------
for i = 1:nr
   
    wh = x_atlas_coverage(i, :) > .5;
    coverage50_index(:, i) = wh';
    
    coverage50_labels{i} = ref_atlas_obj.labels(wh);
    
end


% Make table
% ---------------------------------------------------
Region = atlas_to_parse.labels';
Region_Vol_mm = vol_regions';
Voxels = voxcount_regions';
Ref_region_coverage = round(100 * double(modal_atlas_coverage)); % Percentage of reference atlas regions covered
Atlas_regions_covered = sum(coverage50_index)';

region_table = table(Region, Voxels, Region_Vol_mm, Atlas_regions_covered, modal_label, Ref_region_coverage, modal_atlas_index);

end % function

