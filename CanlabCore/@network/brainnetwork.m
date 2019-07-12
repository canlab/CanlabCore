% A brainnetwork object has these parts:
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
