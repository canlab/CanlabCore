function [region_table, table_legend_text, all_regions_covered, x_counts, x_dice, x_atlas_coverage] = atlas_similarity(atlas_to_parse, ref_atlas_obj)
% Annotate regions in an atlas object with labels from another atlas object
% Take regions in an atlas object (atlas_to_parse) and annotate them with labels and quantitative
% coverage stats from another atlas (ref_atlas_obj)
%
% [region_table, table_legend_text, coverage25_labels, coverage25_index, x_counts, x_dice, x_atlas_coverage] = atlas_similarity(atlas_to_parse, ref_atlas_obj)
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
% all_regions_covered: String matrix list of all reference regions covered by the target region down to 25% coverage, sorted in descending order of coverage
%
% See table_legend_text output for description of table entries.
% to print legend: canlab_print_legend_text(table_legend_text{:})
%
% Tor Wager, July 2018



% Notes:

% Percent of voxels in each atlas region covered by the blob
% Intersection / size of atlas region

% note: vox counts can be zero if regions are (1) outside mask, or (2)
% empty after reslicing to the atlas space

% "coverage" : which atlas regions are covered by the blob? (up to specified % of voxels
% in the atlas region). Large atlas regions/networks will not have high
% coverage unless the blob(s) activate most of the atlas region/network. 
% 
%
% "mode" : what single atlas region best encloses the target blob?
% if blob covers 2 regions completely, the larger is the best match 
%
% "similarity" : jaccard/dice: prioritizes complete coverage of atlas region.
% larger regions will tend to not show up unless the blob
% covers them completely.
%
% could do multi-resolution match with jaccard. small regions would be
% selected if they match, but larger/more general regions would be selected
% if the blob description matches them best.
% e.g., a big region that covers the whole basal ganglia would be labeled
% "BG".  "Mode" would pick the largest subregion. "Coverage" would identify
% multiple subregions. 

%   Notes - from Matlab
%   -----
%   [1]  The Dice similarity coefficient of two sets A and B is
%   expressed as
%
%      dice(A,B) = 2 * |intersection(A,B)| / (|A| + |B|)
%
%   where |A| represents the cardinal of set A. It can also be expressed in
%   terms of true positives (TP), false positives (FP) and false negatives
%   (FN) as
%
%      dice(A,B) = 2 * TP / (2 * TP + FP + FN)
%
%   [2]  The Dice index is related to the Jaccard index according to
%
%      dice(A,B) = 2 * jaccard(A,B) / (1 + jaccard(A,B))
%
%   [1]  The Jaccard similarity coefficient of two sets A and B (also known
%     as intersection over union or IoU) is expressed as
% 
%      jaccard(A,B) = |intersection(A,B)| / |union(A,B)|
% 
%     where |A| represents the cardinal of set A. It can also be expressed in
%     terms of true positives (TP), false positives (FP) and false negatives
%     (FN) as
% 
%        jaccard(A,B) = TP / (TP + FP + FN)

coverage_thresh = .25;      % percent of reference atlas region that must be covered to include it in lists
% Note: if you change this, the help and figure legend and variable names will be wrong/misleading.

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

x_region_coverage = x_counts ./ sum_for_region;

% Modal atlas region for each input region
% ------------------------------------------------------

[max_count_by_region, modal_atlas_index] = max(x_counts, [], 2);

wh_no_atlas_region = max_count_by_region == 0;

modal_atlas_index(wh_no_atlas_region) = 0;  % because max returns 1 where empty counts 

% Modal labels
wh_labels = modal_atlas_index > 0;

modal_label = repmat({'No label'}, nr, 1);
modal_label(wh_labels) = ref_atlas_obj.labels(modal_atlas_index(wh_labels));

modal_label_descriptions = repmat({'No description'}, nr, 1);
modal_label_descriptions(wh_labels) = ref_atlas_obj.label_descriptions(modal_atlas_index(wh_labels))';

if length(modal_label) ~= nr || length(modal_label_descriptions) ~= nr
    warning('atlas_similarity: Label and description lengths do not match number of regions. labels or label_descriptions field in atlas may be incorrectly formatted.');
    error('Quitting.');
end
    
% modal_atlas_coverage: Percent of best reference region covered
%
[modal_atlas_coverage, modal_region_coverage] = deal(zeros(nr, 1));
all_regions_covered = cell(nr, 1);                                  % strings for all regions covered >= 50%

for i = 1:nr
   
    if wh_no_atlas_region(i) 
        % Missing - no atlas regions match
        modal_atlas_coverage(i, 1) = 0;
        modal_region_coverage(i, 1) = 0;
        
    else
        % Coverage of reference atlas region with highest coverage
        mycoverage = x_atlas_coverage(i, :);
        modal_atlas_coverage(i, 1) = mycoverage(modal_atlas_index(i));
        
        % Coverage of blob by best atlas region
        mycoverage = x_region_coverage(i, :);
        modal_region_coverage(i, 1) = mycoverage(modal_atlas_index(i));
        
        % All regions covered in descending order of importance
        [coverage_vals, indx] = sort(x_atlas_coverage(i, :), 'descend');
        wh = coverage_vals >= coverage_thresh;
        
        if sum(wh) == 0
            all_regions_covered{i} = 'No regions';
        end
        
        for j = 1:sum(wh)
           
            all_regions_covered{i} = char(all_regions_covered{i}, sprintf('%s:\t%3.0f%%', ref_atlas_obj.labels{indx(j)}, 100 * coverage_vals(j)));
            
        end

    end
    
end

% Percent in modal atlas region


% Reference regions with at least 25% coverage
% ---------------------------------------------------
for i = 1:nr
   
    wh = x_atlas_coverage(i, :) >= coverage_thresh;
    coverage25_index(:, i) = wh';
    
%     coverage25_labels{i} = ref_atlas_obj.labels(wh);
    
end

% Re-label regions covering many atlas regions
% ---------------------------------------------------
Atlas_regions_covered = sum(coverage25_index)';

wh = Atlas_regions_covered > 5;
modal_label(wh) = {'Multiple regions'};


% Make table
% ---------------------------------------------------

Region = atlas_to_parse.labels';
Region_Vol_mm = vol_regions';
Voxels = voxcount_regions';
Ref_region_perc = round(100 * double(modal_atlas_coverage)); % Percentage of reference atlas regions covered
Perc_covered_by_label = round(100 * double(modal_region_coverage)); % Percentage of reference atlas regions covered

% Atlas_regions_covered = sum(coverage25_index)';

region_table = table(Region, Voxels, Region_Vol_mm, Atlas_regions_covered, modal_label, modal_label_descriptions, Perc_covered_by_label, Ref_region_perc, modal_atlas_index, all_regions_covered);

% Table properties and legend
% ---------------------------------------------------
myvoxsize = sprintf('%d x %d x %d', abs(diag(atlas_to_parse.volInfo.mat(1:3, 1:3))'));


region_table.Properties.Description = sprintf('Regions labeled by reference atlas %s\n', ref_atlas_obj.atlas_name);

region_table.Properties.VariableDescriptions{1} = 'Region: Original region shorttitle';
region_table.Properties.VariableDescriptions{2} = ['Voxels: Number of contiguous ' myvoxsize ' voxels'];
region_table.Properties.VariableDescriptions{3} = 'Region_Vol_mm: Volume of contiguous region in cubic mm.';
region_table.Properties.VariableDescriptions{4} = 'Atlas_regions_covered: Number of reference atlas regions covered at least 25%% by the region. This relates to whether the region covers multiple reference atlas regions';
region_table.Properties.VariableDescriptions{5} = 'Modal_label: Best reference atlas label, defined as reference region with highest number of in-region voxels. Regions covering >25%% of >5 regions labeled as "Multiple regions"';
region_table.Properties.VariableDescriptions{6} = 'Perc_covered_by_label: Percentage of the region covered by the label.';
region_table.Properties.VariableDescriptions{7} = 'Ref_region_perc: Percentage of the label region within the target region.';
region_table.Properties.VariableDescriptions{8} = 'modal_atlas_index: Index number of label region in reference atlas';
region_table.Properties.VariableDescriptions{9} = 'all_regions_covered: All regions covered >5%% in descending order of importance';

% %%%% for compatibility with canlab_print_legend_text
coverage_descrip = {sprintf('For example, if a region is labeled ''TE1a'' and Perc_covered_by_label = 8, Ref_region_perc = 38, and Atlas_regions_covered = 17, this means that 8%%%% of the region''s voxels are labeled TE1a, which is the highest percentage among reference label regions. 38%%%% of the region TE1a is covered by the region. However, the region covers at least 25%%%% of 17 distinct labeled reference regions.\n')};

myrefs = {sprintf('References for atlases:\n')};

myrefs = [myrefs cellstr(unique(ref_atlas_obj.references, 'rows'))'];

% Char array. Now return cells for more flexibility later.
% mystr = unique(ref_atlas_obj.references, 'rows');
% mystr = strvcat(region_table.Properties.Description, mystr);
% mystr = strvcat(mystr, char(region_table.Properties.VariableDescriptions{:}));
% table_legend_text = mystr;

table_legend_text = [region_table.Properties.Description region_table.Properties.VariableDescriptions coverage_descrip myrefs];

% to print: canlab_print_legend_text(table_legend_text{:})

end % function

