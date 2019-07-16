% A brainpathway object has these parts:
% - regions:                A region-class object defining k regions and local patterns in these regions
% - connections:            A series of matrices specifying [k x k] bivariate connections
% - connections.apriori:    [1/0] logical matrix specifying existing connections
% - connections.est:        Estimated connectivity strengths
% - connections.se:         Standard error of estimated connectivity strengths
% - connections.metric      Metric type [r, cos_sim, tau, partial_r]
% - graph_properties: 
% - timeseries_data:        Level-1 (time series, t) cell array with one cell per subject, [t x k] data matrix in each cell, NaN is missing
% - person_data:            Level-2 (person-level, s) data for each region, [s x k] padded matrix, NaN is missing
% - data_properties:        Provenance for what has been done to data

% - processes-metrics


% A brainpathway object has these parts:
% - region_atlas:                An atlas-class object defining k regions
% - WEIGHTS object ->  n fmri_data objects, one per network,
% - connections:            A series of matrices specifying [k x k] bivariate connections
% - connections.apriori:    [1/0] logical matrices specifying existing connections, k x k x n for n networks
% k x n latent variable weights - voxel weights. {1....k} cell k has {v x n}, v voxels x n networks
% OR: weights could be in n fmri_data objects, one per network, with 
% - connections.est:        Estimated connectivity strengths
% - connections.se:         Standard error of estimated connectivity strengths
% - connections.metric      Metric type [r, cos_sim, tau, partial_r]
% - graph_properties: 
% - timeseries_data:        Level-1 (time series, t) cell array with one cell per subject, [t x k] data matrix in each cell, NaN is missing
% - person_data:            Level-2 (person-level, s) data for each region, [s x k] padded matrix, NaN is missing
% - data_properties:        Provenance for what has been done to data

% - processes-metrics

% Extracts local patterns from an image_vector 
% ---------------------------------------------------------
pain_regions_pdm1 = extract_data(pain_regions, pdm1);
% pain_regions_pdm1(1).all_data -> weights are stored in in all_data
% save in .val field, which extract_data will use to extract
for i = 1:length(pain_regions_pdm1)
    pain_regions_pdm1(i).val = pain_regions_pdm1(i).all_data';
    pain_regions_pdm1(i).Z = pain_regions_pdm1(i).all_data;
end
k = length(pain_regions_pdm1);

