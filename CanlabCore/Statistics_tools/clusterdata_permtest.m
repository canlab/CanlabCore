function stats = clusterdata_permtest(data, varargin)
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
%     stats = clusterdata_permtest(data, ['k', # clusters to test, 'nperm,' n_permutations, etc.])
%     stats = clusterdata_permtest(X, 'k', [2:7], 'verbose', false);
%     stats = clusterdata_permtest(X, 'k', [2:7]);
%     stats = clusterdata_permtest(X, 'doplot', false);
%     stats = clusterdata_permtest(X, 'linkage', 'average');
%     stats = clusterdata_permtest(X, 'reducedims', true);
%     stats = clusterdata_permtest(X, 'reducedims', true, 'ndims', 3);
%     stats = clusterdata_permtest(X, 'k', [2:7], 'reducedims', true, 'ndims', 3, 'doplot', false, 'verbose', false);
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
%   **'doplot':** [logical]
%        Flag to control plot output. Default = true.
%
% :Outputs:
%
% :Outputs:
%
%   **stats:** A structure containing the results of the cluster-based permutation test.
%
%       **stats.inputs:** A structure recording the input parameters used in the analysis.
%           - **k:** Vector of candidate numbers of clusters tested.
%           - **distancemetric:** A string specifying the distance metric used.
%           - **linkagemethod:** A string specifying the linkage method for clustering.
%           - **reducedims:** Indicator (flag or value) of whether dimensionality reduction was performed.
%           - **ndims:** Number of dimensions retained after reduction.
%           - **nperm:** Total number of permutations performed.
%           - **verbose:** Logical flag indicating whether verbose output was enabled.
%
%       **stats.cluster_quality_descrip:** A string describing the quality metric used (e.g., 'Mean silhouette value').
%
%       **stats.cluster_quality:** A vector of the mean silhouette values computed for each candidate k.
%
%       **stats.cluster_quality_null_mean:** A vector of the mean silhouette values obtained from the null (permutation) distribution for each candidate k.
%
%       **stats.cluster_quality_null_std:** A vector of the standard deviations of the null silhouette values for each candidate k.
%
%       **stats.cluster_quality_pseudoZ:** A vector of pseudo-Z scores for each candidate k, computed as:
%             (mean silhouette value - null mean) ./ null std.
%
%       **stats.P_val:** A vector of p-values for each candidate k, computed as one minus the proportion of permutations 
%             where the observed mean silhouette value exceeded the null distribution value.
%
%       **stats.max_pseudoZ:** The maximum pseudo-Z value found across candidate k values.
%
%       **stats.wh_best_k:** The index corresponding to the candidate k that produced the maximum pseudo-Z value.
%
%       **stats.best_k:** The candidate number of clusters (k) that yielded the maximum pseudo-Z value.
%
%       **stats.output_table:** A table summarizing the cluster quality metrics across candidate k values.
%           - The table rows are: 'Quality (q)', 'Null q mean', 'Null q std', 'Pseudo-Z', and 'P-value'.
%           - The table columns correspond to the candidate k values (converted to strings).
%           - The table includes a description ('Cluster quality by number of clusters (k)').
%
%       **stats.best_cluster_labels:** The cluster labels (for each observation) corresponding to the best candidate k.
%
%       **stats.all_cluster_labels:** A matrix of cluster labels for all candidate k values.
%
%       **stats.silhouette_values:** The silhouette values for each observation (used in computing cluster quality).
%
%       **stats.linkage_tree:** The linkage tree (hierarchical clustering tree) computed during the clustering analysis.
%
% -------------------------------------------------------------------------
% :Examples:
% ::
%
% rng('default');  % For reproducibility
% n_per_cluster = 33; n_vars = 10;
% X = [gallery('uniformdata',[n_per_cluster n_vars],12); ...
%     gallery('uniformdata',[n_per_cluster n_vars],13)+1.2; ...
%     gallery('uniformdata',[n_per_cluster n_vars],14)+2.5];
% y = [ones(n_per_cluster,1); 2*(ones(n_per_cluster,1)); 3*(ones(n_per_cluster,1))]; % Actual classes
% 
% % Test k = 2 - 7 clusters with 100 random permutations and make a plot of the results: 
% [bestpval,bestmyclass,bestnames,bestX,where,clustnames,stats]=testclustnew(X, [2:7], [], 100);
%
% -------------------------------------------------------------------------
% % Generate null data with no true clusters, and estimate:
% X = randn(n_per_cluster*3, n_vars);
% stats = clusterdata_permtest(X, 'k', [2:7]);
%
% -------------------------------------------------------------------------
% % Load an fmri_data object and cluster the images:
% % Note: you may want to choose the number of dimensions separately.
%
% obj = load_image_set('emotionreg')
% stats1 = clusterdata_permtest(obj.dat', 'k', [2:7], 'reducedims', true);
%
% -------------------------------------------------------------------------
% % Here is a more interesting example of fmri image clustering:
% obj = load_image_set('kragel270')
% obj = rescale(obj, 'zscoreimages'); % images are not on the same scale, so adjust
% plot(obj);
% stats = clusterdata_permtest(obj.dat', 'k', [2:20], 'reducedims', true, 'ndims', 25);
% 
% % Plot the image spatial correlation matrix:
% [cluster_labels_sorted, order_indx] = sort(stats.best_cluster_labels, 'ascend');
% X = X(:, order_indx);
% OUT = plot_correlation_matrix(X, 'doimage', true, 'docircles', false, 'partitions', cluster_labels_sorted);
% % OUT = plot_correlation_matrix(X, 'doimage', true, 'docircles', false, 'partitions', cluster_labels_sorted, 'partitioncolors', seaborn_colors(max(cluster_labels_sorted)));
%
% % Summarize the relationship between the clusters and original study membership
% indic1 = condf2indic(obj.metadata_table.Studynumber);
% indic2 = condf2indic(stats.best_cluster_labels);
% xtab = indic1' * indic2;
% figure; imagesc(xtab);
% xlabel('Clusters'); ylabel('Original studies'); colorbar; 
% colormap(colormap_tor([1 1 1], [1 1 0], [1 .5 0]));
%
% :References:
%   the basic method has been used in these studies (among others):
%   Wager, Scott, & Zubieta 2007, PNAS
%   Wager et al. 2008, Neuron
%   Atlas et al. 2014, Pain
%
% :See also:
%   kmeans, confusionmat, randperm, testclustnew (original function by Chris Summerfield)
%

% -------------------------------------------------------------------------
%     Author and copyright information:
%
%     Copyright (C) 2025  Tor Wager
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

% -------------------------------------------------------------------------
% Defaults and inputs
% -------------------------------------------------------------------------

% defaults not otherwise specified below

k = 2:ceil(size(data, 2)/2);
reducedims = size(data, 2) > size(data, 1) ./ 2; % true if there are <2x as many cases as variables

% Parse input arguments using inputParser and validateattributes.
p = inputParser;
% Required argument: data matrix X.
p.addRequired('data', @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Optional parameters:
p.addParameter('distancemetric', 'euclidean', @(x) validateattributes(x, {'char','string'}, {'nonempty'}));
p.addParameter('linkagemethod', 'ward', @(x) validateattributes(x, {'char','string'}, {'nonempty'}));
p.addParameter('k', k, @(x) validateattributes(x, {'numeric'}, {'vector'}));
p.addParameter('reducedims', reducedims, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('ndims', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));  % empty = choose from data in subfunction
p.addParameter('verbose', true,  @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('doplot', true,  @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('nperm', 100, @(x) validateattributes(x, {'numeric'}, {'scalar'}));

p.parse(data, varargin{:});
ARGS = p.Results;
 
% Distribute parsed parameters back out to variables:
fn = fieldnames(ARGS);

for i = 1:length(fn)
    str = sprintf('%s = ARGS.(''%s'');', fn{i}, fn{i});
    eval(str)
end

% -------------------------------------------------------------------------
% Defaults and inputs
% -------------------------------------------------------------------------

if verbose
    fprintf('Evaluating clustering quality for k values: %s\n', mat2str(k));
    fprintf('Number of permutations: %d\n', nperm);
end

% Cluster with real data
% If we are not permuting, we are done. If we are permuting, repeat with permuted data
% best_cluster_id, best_k, cluster_id are preliminary pre-permutation
% values that will be returned or replaced after permutation
% -------------------------------------------------------------------------
stats = struct;
[stats.best_cluster_labels, ~, all_cluster_labels, mean_silhouette_value, silhouette_values, linkage_tree, ndims] = cluster_dataset(data, varargin{:});

% Keep same ndims for permutation if we are reducing dimensionality
% This avoids problem with permuted data being incommensurate in
% dimensionality with real solution, as permutation disrupts the covariance
% structure. This only matters if reducedims = true, and is ignored
% otherwise

if reducedims

    % add to varargin so we can pass this in; replace default/entered value
    wh_ndims = find(strcmp(varargin, 'ndims'));
    if ~isempty(wh_ndims), varargin(wh_ndims + 1) = {ndims}; end

    if verbose, fprintf('Reduced dimensionality to %3.0f dims\n', ndims); end
end

stats.inputs = struct('k', k, 'distancemetric', distancemetric, 'linkagemethod', linkagemethod, 'reducedims', reducedims, 'ndims', ndims, 'nperm', nperm, 'verbose', verbose);
stats.cluster_quality_descrip = 'Mean silhouette value';

% -------------------------------------------------------------------------
% Cluster with permuted data
% -------------------------------------------------------------------------

numk = numel(k);
null_mean_silhouette_value = zeros(nperm, numk);

if verbose
    % progress bar
    fprintf('Permutations: %03d%%', 0)
end

for i = 1:nperm

    if verbose
        if mod(i, 5) == 0
            fprintf('\b\b\b\b%03d%%', 100 * i ./ nperm)
        end
    end

    data_permuted = permute_X_columns(data);

    % generate a vector of quality scores for each value of k
    [~, ~, ~, null_mean_silhouette_value(i, :)] = cluster_dataset(data_permuted, 'k', k, varargin{:});

end

if verbose
    fprintf(' done! \n')
end

% -------------------------------------------------------------------------
% Save output
% -------------------------------------------------------------------------

nullmean = mean(null_mean_silhouette_value);
nullstd = std(null_mean_silhouette_value);

stats.cluster_quality = mean_silhouette_value;
stats.cluster_quality_null_mean = nullmean;
stats.cluster_quality_null_std = nullstd;
stats.cluster_quality_pseudoZ = (mean_silhouette_value - nullmean) ./ nullstd;
stats.P_val = 1 - sum(mean_silhouette_value > null_mean_silhouette_value) ./ nperm;

[stats.max_pseudoZ, stats.wh_best_k] = max(stats.cluster_quality_pseudoZ);
stats.best_k = k(stats.wh_best_k);

if verbose
    fprintf('Recommended number of clusters: %d (Z-score = %.3f)\n', stats.best_k, stats.max_pseudoZ);
end


% build table
knames = arrayfun(@num2str, stats.inputs.k, 'UniformOutput', false);
t = array2table(mean_silhouette_value(:).');
t.Properties.VariableNames = knames;

t(2, :) = array2table(nullmean(:).');
t(3, :) = array2table(nullstd(:).');
t(4, :) = array2table(stats.cluster_quality_pseudoZ(:).');
t(5, :) = array2table(stats.P_val(:).');
t.Properties.RowNames = {'Quality (q)' 'Null q mean' 'Null q std' 'Pseudo-Z' 'P-value'}';

stats.output_table = t;
t.Properties.Description = 'Cluster quality by number of clusters (k)';

if verbose
    disp(' ')
    disp(t.Properties.Description)
    disp(t)
end

% add additional output
stats.best_cluster_labels = all_cluster_labels(:, stats.wh_best_k);
stats.all_cluster_labels = all_cluster_labels;
stats.silhouette_values = silhouette_values;
stats.linkage_tree = linkage_tree;

if doplot

% Plot a dendrogram for the best clustering solution.
create_figure('cluster_plots', 1, 2); 

dendrogram(gca, stats.linkage_tree, 'ClusterIndices', stats.best_cluster_labels, 'Parent', gca);

subplot(1, 2, 2);

n = size(data, 1);
v = length(unique(stats.best_cluster_labels));  % Number of unique classes

% Run t-SNE to reduce data to two dimensions.
Ytsne = tsne(data);

% Get a color palette with v colors using seaborn_colors from CANlab tools.
colors = seaborn_colors(v);  % palette is a v×3 matrix

% Create a color for each observation based on its label.
pointColors = colors(stats.best_cluster_labels, :);
pointColors = cat(1, pointColors{:});

% Compute marker size proportional to the number of observations.
% (Using an inverse relationship so that larger n yields smaller markers.)
markerSize = 300 / sqrt(n);

% Create the scatter plot.
scatter(Ytsne(:,1), Ytsne(:,2), markerSize, pointColors, 'filled');

% Label the dimensions.
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');

% Set the axis font size.
set(gca, 'FontSize', 18);

% Add a title.
title('t-SNE Visualization');

end


end % main function


% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
% SUBFUNCTIONS
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% PARSE INPUTS
% -------------------------------------------------------------------------

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



% -------------------------------------------------------------------------
% cluster_dataset
% -------------------------------------------------------------------------

function [best_cluster_id, best_k, cluster_id, mean_silhouette_value, silhouette_values, Z, ndims] = cluster_dataset(X, varargin)
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
%        Default = data-dependent (see main function above).
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

verbose = true;

% Parse input arguments using inputParser and validateattributes.
p = inputParser;

k = 2:ceil(size(X, 2)/2);
reducedims = size(X, 2) > size(X, 1) ./ 2; % true if there are <2x as many cases as variables

% Required argument: data matrix X.
p.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Optional parameters:
p.addParameter('distancemetric', 'euclidean', @(x) validateattributes(x, {'char','string'}, {'nonempty'}));
p.addParameter('linkagemethod', 'ward', @(x) validateattributes(x, {'char','string'}, {'nonempty'}));
p.addParameter('k', k, @(x) validateattributes(x, {'numeric'}, {'vector'}));
p.addParameter('reducedims', reducedims, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('ndims', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

% these are not used but keep them because we pass them in from main
% function for convenience
p.addParameter('doplot', false,  @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('verbose', true,  @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('nperm', 100, @(x) validateattributes(x, {'numeric'}, {'scalar'}));

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
        
        try
            ndims = barttest(X, 0.05);
        catch
            ndims = NaN;
        end

        % sometimes this may not work
        % if it doesn't, pick dims with eigenvalues > 1 as a heuristic
        if isnan(ndims)

            [~, ~, latent] = pca(X, 'Economy', true);
            ndims = sum(latent > 1);

        end

        % if verbose, fprintf('Reducing dims: %3.0f\n', ndims); end

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


% -------------------------------------------------------------------------
% permute_X_columns
% -------------------------------------------------------------------------

% Randomly permute each column of the 2-D matrix X independently
function X_permuted = permute_X_columns(X)

[n, k] = size(X);

% Preallocate the output matrix of the same size as X
X_permuted = zeros(n, k);

% Loop through each column
for col = 1:k

    % Randomly permute the rows for the current column
    X_permuted(:,col) = X(randperm(n), col);

end

end % permute

