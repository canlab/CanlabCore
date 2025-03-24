function output = clusterdata_permtest(data, varargin)
% Cluster a dataset and evaluate clustering quality using permutation testing.
%
% This function performs hierarchical clustering on the input dataset X.
% It optionally reduces the dimensionality of the data using PCA if 'reducedims'
% is true (and the number of dimensions is chosen using barttest at p < 0.05 if
% not provided). The function then computes a dendrogram using a specified distance
% metric and linkage method. For a range of cluster numbers (k_values), it computes the
% silhouette values (an index of cluster quality) for each clustering solution and 
% selects the number of clusters (best_k) that maximizes the average silhouette value. 
% The function returns the best cluster assignment, as well as the full set of clustering 
% solutions and quality metrics.
%
% The silhouette values are not uniform across different choices of k, so this function
% uses a permutation-based null hypothesis to test which solution is best. 
% For each k, it computes a clustering quality metric (here defined
% as the sum of within-cluster distances from kmeans) on the actual data and on a
% set of nperm permuted datasets. The separation in quality between the actual
% and permuted solutions is expressed as a Z-score. The function then recommends
% the optimal number of clusters as the value of k with the maximum separation
% (highest Z-score).
%
% :Usage:
% ::
%     output = testclustnew(data, k_values, nperm, [optional inputs])
%
% :Inputs:
%
%   **X:**
%        A numeric matrix [n x p] where each row is an observation and each column is a feature.
%
% :Optional Inputs (passed as name-value pairs):
%
%   **'distancemetric':** [string]
%        The distance metric to use for computing pairwise distances. 
%        Any option available in the pdist.m function is valid.
%        Default = 'euclidean'.
%
%   **'linkagemethod':** [string]
%        The linkage method to use for hierarchical clustering.
%        Any option available in the linkage.m function is valid.
%        Default = 'ward'.
%
%   **'k':** [numeric vector]
%        A vector of cluster numbers to evaluate (e.g., 2:7).
%        Default = 2:7.
%
%   **'reducedims':** [logical]
%        If true, perform dimensionality reduction via PCA before clustering.
%        Default = true if there are <2x as many cases as variables
%
%   **'ndims':** [numeric scalar]
%        The number of dimensions to retain if 'reducedims' is true.
%        If empty, the number of dimensions is chosen using barttest at p < 0.05.
%        Default = [].
%
%   **nperm:**
%        A positive scalar specifying the number of permutations.
%
%   **'verbose':** [logical]
%        Flag to control verbose output. Default = true.
%
% :Outputs:
%
%   **best_cluster_id:**
%        A numeric vector containing the cluster assignment (for best_k) for each observation.
%
%   **best_k:**
%        The number of clusters with the maximum average silhouette value.
%
%   **cluster_id:**
%        A matrix where each column corresponds to the cluster assignments for each value in k.
%
%   **mean_silhouette_value:**
%        A vector of the mean silhouette values for each clustering solution in k.
%
%   **silhouette_values:**
%        A matrix of silhouette values for each observation (rows) and each k (columns).
%
% :Outputs:
%
%   **output:**
%        A structure containing:
%           - k_values: the tested k values.
%           - actual_quality: vector of clustering quality metrics for the actual data.
%           - null_quality: matrix of clustering quality metrics from permutations
%                           (size: length(k_values) x nperm).
%           - z_scores: vector of Z-scores representing the separation between
%                       actual and null distributions.
%           - recommended_k: the value of k with the maximum Z-score.
%
% :Examples:
% ::
%     data = rand(100, 5);
%     k_vals = 2:10;
%     nperm = 100;
%     results = testclustnew(data, k_vals, nperm, 'verbose', true, 'distance', 'sqEuclidean');
%
% :References:
%   (Add any relevant references here.)
%
% :See also:
%   kmeans, confusionmat, randperm, testclustnew (original function by Chris Summerfield)
%
% -------------------------------------------------------------------------
%     Author and copyright information:
%
%     Copyright (C) 2025  Your Name
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
% -------------------------------------------------------------------------

% defaults
k_values = 2:ceil(size(data, 2)/2);
nperm = 100;
reducedims = size(data, 2) > size(data, 1) ./ 2; % true if there are <2x as many cases as variables
verbose = true;


% Parse input arguments using inputParser and validateattributes.
p = inputParser;
% Required argument: data matrix X.
p.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Optional parameters:
p.addParameter('distancemetric', 'euclidean', @(x) validateattributes(x, {'char','string'}, {'nonempty'}));
p.addParameter('linkagemethod', 'ward', @(x) validateattributes(x, {'char','string'}, {'nonempty'}));
p.addParameter('k', k_values, @(x) validateattributes(x, {'numeric'}, {'vector'}));
p.addParameter('reducedims', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('ndims', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('verbose', true,  {'logical'}, {'scalar'});

p.parse(X, varargin{:});
ARGS = p.Results;
 
% Distribute parsed parameters.
distancemetric = ARGS.distancemetric;
linkagemethod = ARGS.linkagemethod;
k = ARGS.k;
reducedims = ARGS.reducedims;
ndims = ARGS.ndims;

if verbose
    fprintf('Evaluating clustering quality for k values: %s\n', mat2str(k_values));
    fprintf('Number of permutations: %d\n', nperm);
end



% Plot a dendrogram for the best clustering solution.
figure; 
dendrogram(Z, 'ClusterIndices', best_cluster_id);






numK = numel(k_values);
actual_quality = zeros(numK, 1);
null_quality = zeros(numK, nperm);

% Loop over each k value.
for i = 1:numK
    k = k_values(i);
    % Perform clustering using kmeans. Replicates ensure stability.
    [~, ~, sumd] = kmeans(data, k, 'Distance', distanceMetric, 'Replicates', 5);
    actual_quality(i) = sum(sumd);
    
    % Permutation testing: shuffle the rows of data.
    for p = 1:nperm
        permutedData = data(randperm(size(data,1)), :);
        [~, ~, sumd_perm] = kmeans(permutedData, k, 'Distance', distanceMetric, 'Replicates', 5);
        null_quality(i, p) = sum(sumd_perm);
    end
    
    if verbose
        fprintf('k = %d: Actual quality = %.3f\n', k, actual_quality(i));
    end
end

% Compute Z-scores for each k: (actual - mean(null))/std(null)
z_scores = zeros(numK, 1);
for i = 1:numK
    mu_null = mean(null_quality(i, :));
    sigma_null = std(null_quality(i, :));
    z_scores(i) = (actual_quality(i) - mu_null) / sigma_null;
end

% Determine the recommended k: maximum Z-score.
[~, bestIdx] = max(z_scores);
recommended_k = k_values(bestIdx);

% Create output structure.
output = struct;
output.k_values = k_values;
output.actual_quality = actual_quality;
output.null_quality = null_quality;
output.z_scores = z_scores;
output.recommended_k = recommended_k;

if verbose
    fprintf('Recommended number of clusters: %d (Z-score = %.3f)\n', recommended_k, z_scores(bestIdx));
end

end

function ARGS = parse_inputs(varargin)
    p = inputParser;
    % Required inputs
    % p.addRequired('data', @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
    p.addRequired('k_values', @(x) validateattributes(x, {'numeric'}, {'vector', '>=', 2}));
    p.addRequired('nperm', @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    
    % Optional parameters
    p.addParameter('verbose', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    p.addParameter('distance', 'sqEuclidean', @(x) ischar(x));
    
    p.parse(varargin{:});
    ARGS = p.Results;
end




function [best_cluster_id, best_k, cluster_id, mean_silhouette_value, silhouette_values, Z] = cluster_dataset(X, varargin)
% cluster_dataset Perform hierarchical clustering on dataset X.
%
% This function performs hierarchical clustering on the input dataset X.
% It optionally reduces the dimensionality of the data using PCA if 'reducedims'
% is true (and the number of dimensions is chosen using barttest at p < 0.05 if
% not provided). The function then computes a dendrogram using a specified distance
% metric and linkage method. For a range of cluster numbers (k), it computes the
% silhouette values for each clustering solution and selects the number of clusters
% (best_k) that maximizes the average silhouette value. The function returns the best
% cluster assignment, as well as the full set of clustering solutions and quality metrics.
%
% :Usage:
% ::
%     [best_cluster_id, best_k, cluster_id, mean_silhouette_value, silhouette_values] = ...
%             cluster_dataset(X, 'distancemetric', 'euclidean', 'linkagemethod', 'ward', ...
%                             'k', 2:7, 'reducedims', true, 'ndims', []);
%
% :Inputs:
%
%   **X:**
%        A numeric matrix [n x p] where each row is an observation and each column is a feature.
%
% :Optional Inputs (passed as name-value pairs):
%
%   **'distancemetric':** [string]
%        The distance metric to use for computing pairwise distances.
%        Default = 'euclidean'.
%
%   **'linkagemethod':** [string]
%        The linkage method to use for hierarchical clustering.
%        Default = 'ward'.
%
%   **'k':** [numeric vector]
%        A vector of cluster numbers to evaluate (e.g., 2:7).
%        Default = 2:7.
%
%   **'reducedims':** [logical]
%        If true, perform dimensionality reduction via PCA before clustering.
%        Default = true.
%
%   **'ndims':** [numeric scalar]
%        The number of dimensions to retain if 'reducedims' is true.
%        If empty, the number of dimensions is chosen using barttest at p < 0.05.
%        Default = [].
%
% :Outputs:
%
%   **best_cluster_id:**
%        A numeric vector containing the cluster assignment (for best_k) for each observation.
%
%   **best_k:**
%        The number of clusters with the maximum average silhouette value.
%
%   **cluster_id:**
%        A matrix where each column corresponds to the cluster assignments for each value in k.
%
%   **mean_silhouette_value:**
%        A vector of the mean silhouette values for each clustering solution in k.
%
%   **silhouette_values:**
%        A matrix of silhouette values for each observation (rows) and each k (columns).
%
% :Examples:
% ::
%     X = randn(100, 20); % Generate a random dataset with 100 observations and 20 features.
%     [best_cluster_id, best_k, cluster_id, mean_silhouette_value, silhouette_values] = ...
%         cluster_dataset(X, 'distancemetric', 'euclidean', 'linkagemethod', 'ward', ...
%                         'k', 2:7, 'reducedims', true, 'ndims', []);

% Parse input arguments using inputParser and validateattributes.
p = inputParser;
% Required argument: data matrix X.
p.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Optional parameters:
p.addParameter('distancemetric', 'euclidean', @(x) validateattributes(x, {'char','string'}, {'nonempty'}));
p.addParameter('linkagemethod', 'ward', @(x) validateattributes(x, {'char','string'}, {'nonempty'}));
p.addParameter('k', 2:7, @(x) validateattributes(x, {'numeric'}, {'vector'}));
p.addParameter('reducedims', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('ndims', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.parse(X, varargin{:});
ARGS = p.Results;
 
% Distribute parsed parameters.
distancemetric = ARGS.distancemetric;
linkagemethod = ARGS.linkagemethod;
k = ARGS.k;
reducedims = ARGS.reducedims;
ndims = ARGS.ndims;

% If dimensionality reduction is requested, reduce dimensions using PCA.
if reducedims
    % If ndims is empty, determine the number of dimensions using barttest at p < 0.05.
    if isempty(ndims)
        ndims = barttest(X, 0.05);
    end
    % Perform PCA retaining the specified number of dimensions.
    [~, X_to_cluster] = pca(X, 'Economy', true, 'NumComponents', ndims);
else
    X_to_cluster = X;
end

% Compute the pairwise distances using the specified distance metric.
Y = pdist(X_to_cluster, distancemetric);
% Perform hierarchical clustering using the specified linkage method.
Z = linkage(Y, linkagemethod);

% Obtain clustering solutions for each value in k.
cluster_id = cluster(Z, "maxclust", k);

% Compute the silhouette values for each clustering solution.
for i = 1:length(k)
    silhouette_values(:, i) = silhouette(X_to_cluster, cluster_id(:, i), distancemetric);
end

% Compute the mean silhouette value for each k.
mean_silhouette_value = mean(silhouette_values);

% Select the best k as the one with the maximum mean silhouette value.
[~, wh_best] = max(mean_silhouette_value);
best_k = k(wh_best);
best_cluster_id = cluster_id(:, wh_best);


end  % cluster_dataset
