tic

atlas_obj = load_atlas('canlab2018_2mm', 'noverbose');

%%
r = region2atlas(b_regions_fdr);

%k_orig = num_regions(r); % original regions

%r = resample_space(r, atlas_obj);

% Resample atlas space so they match
% Will compress indices and remove regions that do not match
%atlas_obj = resample_space(atlas_obj, r); % resample atlas to region space - really slow if atlas_obj is higher-res

[region_table, coverage50_labels, coverage50_index, x_counts, x_dice, x_atlas_coverage] = atlas_similarity(r, atlas_obj);

toc



% if k ~= k_orig
%     warning('Some regions lost in resampling to atlas!');
% end


%k = length(a_regions_fdr);


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


