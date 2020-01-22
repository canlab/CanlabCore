function fmri_dat_obj = brainpathway2fmri_data(obj, data_to_add)
% Extract data from any field of a brainpathway object (or external data
% vector) and save it as an fmri_data object for viewing and statistics.
%
% fmri_dat_obj = brainpathway2fmri_data(bs, data matrix or vector)
%
% - Uses obj.region_atlas to define regions
% - Data should either be (v) voxels x images or (k) regions x images
% - If data is region-level (k-length), will be expanded to voxel-level in images
% - 
% Examples:
% ------------------------------------------------------------------
% Load brainpathway_data_gsr_censoring_pipeline.mat % from OLP/Yoni
% bs = bct_toolbox_undirected_graph_metrics(bs);
% mean_degree = nanmean(bs.graph_properties.regions.degree')'; % mean degree for each node
%
% % Plot montage of mean degree on brain
% fmri_dat_obj = brainpathway2fmri_data(bs, bs.graph_properties.regions.degree);
% montage(mean(fmri_dat_obj));

v = obj.region_atlas.volInfo.n_inmask;  % voxels
k = num_regions(obj.region_atlas);      % regions

switch size(data_to_add, 1)
    case v
        % do nothing, OK
    case k
        
        % Need to expand to correlation value for each voxel
        data_to_add = expand_values_region2voxel(obj, data_to_add);
        
    otherwise
        error('brainpathway2fmri_data: size(data_to_add, 1) must = num voxels or regions in obj.region_atlas');
end

fmri_dat_obj = fmri_data(image_vector('volInfo', obj.region_atlas.volInfo, 'dat', data_to_add));

%fmri_dat_obj.image_names = char(obj.region_atlas.labels{to_extract});
fmri_dat_obj.source_notes = 'Maps extracted with brainpathway2fmri_data';

end % function


