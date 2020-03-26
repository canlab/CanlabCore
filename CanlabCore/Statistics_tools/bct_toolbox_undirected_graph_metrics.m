function [graph_prop, graph_prop_glob] = bct_toolbox_undirected_graph_metrics(r, varargin)
% Calculate some Sporns BCT toolbox functions
%
% graph_prop = bct_toolbox_undirected_graph_metrics(r, [threshold_input])
% 
% Inputs:
% r = correlation matrix
% threshold_input = logical positive sig matrix or %linkdensity (scalar).
%       Default = 10 (top 10% of links retained)
%       If sig matrix includes sig negative associations, results could be
%       misleading; they will be considered to be connected just as
%       positive associations are. Best to use only sig positive
%       associations.
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
% 
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

% Prep r for BCT undirected 
r = double(r);
r = (r' + r) ./ 2;          % enforce symmetry (rounding error possible)
r = r - eye(size(r));       % for BCT and squareform

% Threshold: Use sig matrix or link density
sigmat = [];
if isempty(varargin)
    thr_percent_dens = 10; % 10% link density
elseif ismatrix(varargin{1})
    thr_percent_dens = NaN;
    sigmat = varargin{1};
else
    thr_percent_dens = varargin{1};
end

if isempty(sigmat) 
    % use thr for link density
    thr_val = prctile(squareform(r), 100 - thr_percent_dens); % uses stats toolbox
    sigmat = r > thr_val;                                     % positive connections only
end

% validate sigmat: logical
validateattributes(sigmat,{'logical'}, {'nonempty' 'nonnegative'});

bu_matrix = sigmat; 
wu_matrix = r .* sigmat;

% Node-level properties: Undirected binary and weighted networks
graph_prop = table();

graph_prop.community = modularity_und(bu_matrix, 1);           % C: communities from adjacency matrix
graph_prop.community_w = modularity_und(wu_matrix, 1);
graph_prop.community_w2 = community_louvain(wu_matrix, 1, [], 'negative_asym');
graph_prop.core_w = core_periphery_dir(wu_matrix, 1)';         % Core-periphery

graph_prop.degree = degrees_und(bu_matrix)';                   % Node degree
graph_prop.strength = strengths_und(wu_matrix)';               % weighted strength. tested: same results as using r_mean
graph_prop.betweenness = betweenness_bin(bu_matrix)';
graph_prop.clustercoef = clustering_coef_bu(bu_matrix);        % Clustering coefficient
graph_prop.local_efficiency = efficiency_bin(bu_matrix, 1);    % Local efficiency

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
