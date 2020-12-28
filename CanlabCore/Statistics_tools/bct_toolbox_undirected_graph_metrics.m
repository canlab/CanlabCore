function [graph_prop, graph_prop_glob] = bct_toolbox_undirected_graph_metrics(r, thresh, varargin)
% Calculate some Sporns BCT toolbox functions
%
% graph_prop = bct_toolbox_undirected_graph_metrics(r, [threshold_input])
% 
% Inputs:
% r = correlation matrix
% thresh = 0 to 1 value. (.1 is a common value)
%
% optional: 
%   'doplots' - show matrices
%   'doweighted' -- this increase compute time by several orders of
%   magnitude. Default = false
%
% Outputs:
% graph_prop = A table of node-level graph metrics 
% graph_prop_glob = A table of global graph metrics
% 
% For BCT toolbox, see:
% https://sites.google.com/site/bctnet/
%
% For descriptions of metrics, see:
% https://sites.google.com/site/bctnet/measures/list
%
% Notes:
% - Some metrics are calculated on thresholded, weighted (continuous valued) matrix r
%   Sig matrix or link density threshold will be applied. 
%   To use without any thresholding, sig mat should be all ones, or link density 100
% - Others are calculated on thresholded, binarized matrix (positive connections only)
%   If you enter a sig matrix, this will be used as the input to graph metric functions
%   If you enter a link density, a binary matrix will be calculated
% - Fisher r to Z is computed, which will have little impact in many cases
%   but may help in some cases
% - Many BCT functions will not use negative values and ignore weights, so require thresholding. 
%   P-values/statistical significance is one way to threshold.
%   Another way would be using multiple arbitrary thresholds (e.g., 10% link density)
%   The latter may be practical for individual subjects, unless individual-subject stats are calculated and saved
% - This function prioritizes 
% - Degree, etc. may possibly be related to confounds (head motion), and may want to adjust/test for this.
%
% Examples:
% -------------------------------------------------------
% % start with r, a series of correlation matrices, one per subject (k x k x n_subjects)
% OUT = ttest3d(r);                                             % Get some basic stats
% graph_prop = bct_toolbox_undirected_graph_metrics(OUT.r);     % input group mean correlation matrix
%
% rr = corr(table2array(graph_prop));                           % correlate the metrics
% plot_correlation_matrix(rr, 'names', graph_prop.Properties.VariableNames)
%
% % Note: correlations among community vectors are not meaningful. 
% % for similarity among communities/modules detected with different methods, see partition_distance.m

% Check path
if isempty(which('threshold_proportional.m'))
    error('The BCT toolbox does not seem to be on your Matlab path. See https://sites.google.com/site/bctnet/');
end

doplots = false;
doweighted = false;
if any(strcmp(varargin, 'doplots')), doplots = true; end
if any(strcmp(varargin, 'doweighted')), doweighted = true; end


% Prep r for BCT undirected 
r = double(r);
r = (r' + r) ./ 2;          % enforce symmetry (rounding error possible)
r = r - eye(size(r));       % for BCT and squareform

% Missing regions/constant values will give NaNs, so need to account for them

numnans = any(isnan(r));
if sum(numnans) > .10 * length(numnans)
    warning('More than 10% of vars have NaN correlation values [missing data?]. Replacing NaNs with 0s. Be careful about effects on subsequent metrics');
end

r(isnan(r)) = 0;

% Threshold: Use sig matrix or link density
bu_matrix = weight_conversion(threshold_proportional(r, thresh), 'binarize');

% Make weighted matrix based on Fisher transform
z = rToZ(r);
z(isinf(z)) = max(z(~isinf(z))); % correlations of 1.00 get transformed to Inf, which causes BCT to crash in some functions. Replace w/ max value. Note that this code does not correctly handle case of r = -1.00, as this is very unlikely.
wu_matrix = z;

% view
if doplots
    create_figure('BCT networks', 1, 2)
    imagesc(bu_matrix), colorbar
    subplot(1,2,2)
    imagesc(wu_matrix), colorbar
end

% Node-level properties: Undirected binary and weighted networks
graph_prop = table();

graph_prop.community = modularity_und(bu_matrix, 1);           % C: communities from adjacency matrix
graph_prop.core_w = core_periphery_dir(wu_matrix, 1)';         % Core-periphery

graph_prop.degree = degrees_und(bu_matrix)';                   % Node degree
graph_prop.betweenness_bin = betweenness_bin(double(bu_matrix))';  % note: logical did not work in some cases...
graph_prop.clustercoef_bin = clustering_coef_bu(bu_matrix);        % Clustering coefficient
graph_prop.local_efficiency_bin = efficiency_bin(bu_matrix, 1);    % Local efficiency

% weighted matrix computations
if doweighted
    graph_prop.strength = strengths_und(wu_matrix)';               % weighted strength. tested: same results as using r_mean
    graph_prop.local_efficiency_weighted = efficiency_wei(wu_matrix, 1);    % Local efficiency
    graph_prop.eigenvector_centrality = eigenvector_centrality_und(wu_matrix);    % Eigenvector centrality (similar to PageRank)
    graph_prop.clustercoef_weighted = clustering_coef_wu(wu_matrix);        % Clustering coefficient
    graph_prop.betweenness_weighted = betweenness_wei(wu_matrix);  % note: slow
    
    graph_prop.community_w = modularity_und(wu_matrix, 1);
    graph_prop.community_w2 = community_louvain(wu_matrix, 1, [], 'negative_asym');
end
    % Others to consider
% assortativity_bin(bu_network, 0);
% kcore_bu
% graph_prop.betweenness_w = betweenness_wei(OUT.r);            % slow
% graph_prop.rich_club_coeff = rich_club_bu(bu_network)';       % output len not same as number of nodes

% Note: for similarity among communities/modules, see partition_distance.m

% Global properties
graph_prop_glob = table();

graph_prop_glob.trans_bu = transitivity_bu(bu_matrix);
graph_prop_glob.glob_efficiency = efficiency_bin(bu_matrix);

end
