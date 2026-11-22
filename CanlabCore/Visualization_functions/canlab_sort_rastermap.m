function [isort, cluster_labels, OUT] = canlab_sort_rastermap(X, varargin)
% Sort the rows of a data matrix so that rows with similar activity are adjacent (Rastermap; Stringer et al. 2024)
%
% canlab_sort_rastermap is a MATLAB implementation of the Rastermap algorithm
% (Stringer, Zhong, Syeda, Du, Kesa & Pachitariu, 2024, Nature Neuroscience),
% which orders the rows of a [rows x columns] matrix -- e.g., neurons x time,
% voxels x images, or voxels x time points -- along a 1-D embedding so that rows
% with similar activity patterns sit next to each other. Displaying the matrix
% in this order (a "rastermap" or sorted carpet plot) reveals structure that is
% invisible in the original row order: groups of co-active rows, graded
% transitions between them, and (with the 'sequence' preset) sequences of
% activity that travel across rows over time.
%
% The algorithm follows the reference Python implementation
% (https://github.com/MouseLand/rastermap) step by step:
%
%   1. Normalize: z-score each row across columns, then project out the mean
%      across rows at each column (so a global signal does not dominate).
%   2. Reduce: keep the top n_PCs principal components of the rows.
%   3. Cluster: "scaled" k-means (correlation distance, k-means++ seeding) of the
%      rows into n_clusters clusters in PC space.
%   4. Sort clusters: order the clusters to maximize the match between their
%      similarity matrix and a template matrix that rewards placing similar
%      clusters close together. The template combines a global term
%      (-log of the distance between positions in the ordering) and, with
%      'locality' > 0, a local term (only immediate neighbours count). The
%      search repeatedly moves (and optionally reverses) segments of the
%      ordering, keeping the best move, until no move improves the match.
%   5. Upsample: interpolate the sorted cluster centers onto a finer grid
%      (grid_upsample x n_clusters nodes) and assign each row to the grid node
%      it correlates with best. The grid position is the row's 1-D embedding,
%      and isort sorts rows by it.
%
% This port computes the score of every candidate segment move in closed form
% with batched matrix products instead of re-scoring each permuted matrix, so
% the cluster-sorting step runs quickly in plain MATLAB (no MEX/compiled code).
% Results match the reference implementation up to the random k-means
% initialization (see Programmers' notes).
%
% Optionally, canlab_sort_rastermap also groups the rows into k hierarchical clusters
% (Ward linkage on the upsampled grid nodes) and returns a cluster label for
% every row, and it can plot the sorted matrix with clusters in different
% colours.
%
% :Usage:
% ::
%
%     [isort, cluster_labels, OUT] = canlab_sort_rastermap(X, [optional inputs])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2026 Tor Wager
%     MATLAB implementation of Rastermap, Copyright (C) 2023 Howard Hughes
%     Medical Institute, authored by Carsen Stringer and Marius Pachitariu
%     (GPL-3.0; https://github.com/MouseLand/rastermap).
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
% ..
%
% :Inputs:
%
%   **X:**
%        [n_rows x n_columns] numeric matrix. Rows are the things to be sorted
%        (neurons, voxels, regions); columns are the samples they are compared
%        across (time points, images, trials, conditions). For an fmri_data
%        object obj, use X = obj.dat to sort voxels by their profile across
%        images. Rows containing NaN or with zero variance are not sorted;
%        they are placed at the end of isort.
%
% :Optional Inputs:
%
%   **'n_clusters', [integer]:**
%        Number of k-means clusters that are sorted and then upsampled;
%        default = 100. Must be < 200. Values above ~150 make the sort slow.
%        Data with <= 50 valid rows skip clustering and sort the rows directly.
%
%   **'n_PCs', [integer]:**
%        Number of principal components used; default = 200 (capped at the
%        number of rows/columns).
%
%   **'locality', [0 to 1]:**
%        How local the sort is; default = 0. 0 = global sort (the overall
%        arrangement of clusters matters); 1 = highly local (only immediate
%        neighbours matter), good for finding sequences.
%
%   **'time_lag_window', [integer]:**
%        Number of columns (time points) into the future over which cluster
%        cross-correlations are computed (the maximum over lags is used);
%        default = 0. Use > 0 for sequence finding in time series.
%
%   **'symmetric', [logical]:**
%        Use both positive and negative time lags; default = false. Only
%        matters when time_lag_window > 0.
%
%   **'sequence':**
%        Keyword. Preset for data with sequences: time_lag_window = 10,
%        locality = 0.75 (as Rastermap's sequence_settings).
%
%   **'grid_upsample', [integer]:**
%        Upsampling factor for the cluster grid; default = 10. 0 = no
%        upsampling (rows are ordered by their k-means cluster only).
%
%   **'time_bin', [integer]:**
%        Average this many adjacent columns before the PCA; default = 0
%        (no binning).
%
%   **'normalize', [logical]:**
%        z-score each row before sorting; default = true.
%
%   **'mean_time', [logical]:**
%        Project out the across-row mean at each column; default = true.
%
%   **'bin_size', [integer]:**
%        Number of sorted rows averaged into each row of OUT.X_embedding (the
%        display matrix); default = 0, which picks max(1, floor(n_rows/500)).
%
%   **'n_iter', [integer]:**
%        Maximum number of moves in each pass of the cluster sort;
%        default = 400.
%
%   **'random_state', [integer]:**
%        Seed for the random number generator (k-means seeding); default = 0.
%        The global RNG state is restored afterwards.
%
%   **'kmeans_init_max', [integer]:**
%        Maximum number of rows used to choose the k-means++ seeds;
%        default = 20000. Seeding cost grows with n_rows x n_clusters x 100,
%        so a random subset keeps large datasets (e.g., whole-brain voxels)
%        fast. All rows are still used in the k-means iterations.
%
%   **'kmeans_max_samples', [integer]:**
%        Maximum number of rows used for the k-means iterations; default = Inf
%        (all rows, as in the reference implementation). Set to e.g. 50000
%        for interactive use on very large matrices: k-means then runs on a
%        random subset and every row is assigned afterwards.
%
%   **'k', [integer]:**
%        Number of row clusters to return in cluster_labels (and to colour
%        in the plot); default = 0 (no clustering). Also accepted as
%        'cluster', k.
%
%   **'cluster_method', ['hierarchical' or 'contiguous']:**
%        How rows are grouped when k > 0; default = 'hierarchical'.
%        'hierarchical': Ward linkage on the upsampled grid nodes (each row
%        inherits the label of its grid node). Clusters can be split across
%        the sorted order if similar patterns appear in several places.
%        'contiguous': Ward-style agglomeration restricted to neighbouring
%        positions in the rastermap order, so each cluster is one contiguous
%        band of rows in the sorted display.
%
%   **'colors':**
%        Colours for plotting clusters: an [n x 3] RGB matrix, a colormap
%        name (e.g., 'lines', 'turbo', 'parula'), or a cell array of [r g b]
%        triplets. Default = seaborn_colors(k).
%
%   **'clim', [low high]:**
%        Colour limits for the plotted (z-scored) X_embedding; default =
%        [0 1.5], as in the Rastermap paper and GUI.
%
%   **'plot' / 'doplot', [logical]:**
%        Plot the sorted matrix; default = false. Keyword 'plot' turns it on.
%        With k > 0, each cluster is drawn in its own colour, with a colour
%        strip marking cluster membership along the rows.
%
%   **'verbose' / 'doverbose', [logical]:**
%        Print progress and timing; default = true. Keyword 'noverbose' turns
%        it off.
%
% :Outputs:
%
%   **isort:**
%        [n_rows x 1] row order. X(isort, :) is the sorted matrix. Invalid
%        rows (NaN or zero variance) are placed at the end.
%
%   **cluster_labels:**
%        [n_rows x 1] cluster label (1..k) for each row of X (in the original
%        row order), numbered in order of first appearance from the top of
%        the sorted display. NaN for invalid rows. Empty if k = 0.
%
%   **OUT:**
%        Structure with fields:
%          embedding       - [n_rows x 1] 1-D position of each row (NaN if invalid)
%          isort           - row order (same as first output)
%          igood           - logical [n_rows x 1], rows that were sorted
%          cluster_labels  - same as second output
%          embedding_clust - [n_good x 1] sorted k-means cluster of each valid row
%          X_embedding     - [n_bins x n_columns] sorted, row-binned, z-scored
%                            data, for display
%          bin_size        - rows per bin in X_embedding
%          bin_labels      - cluster label of each X_embedding row (if k > 0)
%          U_nodes         - [n_clusters x n_PCs] sorted cluster centers (PC space)
%          X_nodes         - [n_clusters x n_columns] cluster centers as profiles
%          cc              - [n_clusters x n_clusters] sorted cluster similarity
%          U_upsampled     - [n_grid x n_PCs] upsampled grid nodes
%          linkage_tree    - Ward linkage tree (k > 0, 'hierarchical')
%          settings        - all settings used
%          timing          - seconds spent in each step, and total
%          fig_handle      - figure handle, if plotted
%
% :Examples:
% ::
%
%    % 1. Synthetic data: 2000 rows tuned to positions on a ring, shuffled.
%    %    The rastermap order recovers the ring; hierarchical order does not
%    %    preserve the gradual transitions as well.
%    n = 2000; t = 600;
%    phase = 2*pi*rand(n, 1);                       % latent position of each row
%    stim = 2*pi*cumsum(randn(1, t))/20;            % slowly drifting stimulus
%    X = exp(2*cos(phase - stim)) + randn(n, t);    % tuned responses + noise
%    [isort, ~, OUT] = canlab_sort_rastermap(X, 'plot');
%    figure; plot(phase(isort), '.'); ylabel('true phase'); xlabel('sorted row')
%
%    % 2. Sort voxels of an fMRI dataset by their profile across images and
%    %    group them into 6 clusters, plotted in different colours.
%    obj = load_image_set('emotionreg', 'noverbose');
%    [isort, labels] = canlab_sort_rastermap(obj.dat, 'k', 6, 'plot');
%
%    % 3. Same, with contiguous bands and a custom colormap
%    [isort, labels] = canlab_sort_rastermap(obj.dat, 'k', 6, 'cluster_method', 'contiguous', ...
%        'colors', 'turbo', 'plot');
%
%    % 4. Map the cluster labels back onto the brain
%    lab_img = mean(obj); lab_img.dat = labels; montage(lab_img, 'indexmap', seaborn_colors(6))
%
%    % 5. Sequences in time series data (e.g., neurons x time): use the preset
%    [isort, ~, OUT] = canlab_sort_rastermap(X, 'sequence', 'plot');
%
%    % 6. Sorted carpet plot in fmri_data.plot
%    plot(obj, 'rastermap')
%
% :References:
%   Stringer C, Zhong L, Syeda A, Du F, Kesa M, Pachitariu M (2024).
%   Rastermap: a discovery method for neural population recordings.
%   Nature Neuroscience. https://doi.org/10.1038/s41593-024-01783-4
%
%   Reference implementation: https://github.com/MouseLand/rastermap
%
%   Arthur D, Vassilvitskii S (2007). k-means++: the advantages of careful
%   seeding. ACM-SIAM Symposium on Discrete Algorithms.
%
% :See also:
%   - canlab_sort_distance_matrix (hierarchical-clustering order of a distance matrix)
%   - plot_correlation_matrix ('reorder_by_clustering'), clusterdata_permtest
%   - fmri_data.plot ('rastermap' option), seaborn_colors
%

% ..
%    Programmers' notes:
%    Nov 2026  Tor Wager (with Claude Code): MATLAB port of rastermap v1.0
%    (rastermap.py, cluster.py, sort.py, upsample.py). Differences from the
%    reference:
%      - The traveling-salesman cluster sort scores all segment moves for a
%        given segment length in closed form (score = within-remainder terms
%        via 2-D cumulative sums + segment/remainder cross terms via paged
%        matrix products and diagonal sums) rather than permuting and
%        re-scoring each candidate matrix. The search order, move set,
%        acceptance threshold (1e-3), coarse (n_skip) and fine passes, and
%        tie-breaking follow tsp_fast() exactly, so for a given similarity
%        matrix the sort is the same up to floating point.
%      - Random numbers come from MATLAB's RNG, so k-means seeds differ from
%        numpy's; final orders differ accordingly (as with any seed change).
%      - k-means++ seeding can use a random subset of rows (kmeans_init_max)
%        and k-means can run on a subset (kmeans_max_samples) for speed.
%      - The PCA is exact (eigendecomposition of the smaller Gram matrix)
%        rather than randomized TruncatedSVD.
%      - n_splits (split-and-resort refinement) and user-supplied BBt / U_nodes
%        are not implemented.
%      - Hierarchical clustering of rows, and plotting, are CANlab additions.
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

% Parse special command keywords and remove them before inputParser

doplot = false;
plot_idx = strcmpi(varargin, 'plot');
if any(plot_idx)               % Override: omit 'doplot' key/value pair
    doplot = true;
    varargin(plot_idx) = [];   % remove so inputParser doesn't see it
end

verbose = true;
verbose_idx = strcmpi(varargin, 'noverbose');
if any(verbose_idx)
    verbose = false;
    varargin(verbose_idx) = [];   % remove so inputParser doesn't see it
end

% 'sequence' preset (Rastermap's sequence_settings): sets defaults that
% explicit key/value pairs below can still override
locality_default = 0;
time_lag_default = 0;
seq_idx = strcmpi(varargin, 'sequence');
if any(seq_idx)
    locality_default = 0.75;
    time_lag_default = 10;
    varargin(seq_idx) = [];
end

% 'cluster', k is an alias for 'k', k
wc = find(strcmpi(varargin, 'cluster'));
for i = wc(:)', varargin{i} = 'k'; end

% Use inputParser to parse key/value pairs
% First add obligatory/non-conditional keywords

valfcn_posint  = @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
valfcn_logical = @(x) islogical(x) || isnumeric(x);

p = inputParser;
p.addRequired('X', @(x) validateattributes(x, {'numeric', 'logical'}, {'2d', 'nonempty'}));
p.addParameter('n_clusters', 100, valfcn_posint);
p.addParameter('n_PCs', 200, valfcn_posint);
p.addParameter('locality', locality_default, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
p.addParameter('time_lag_window', time_lag_default, valfcn_posint);
p.addParameter('symmetric', false, valfcn_logical);
p.addParameter('grid_upsample', 10, valfcn_posint);
p.addParameter('time_bin', 0, valfcn_posint);
p.addParameter('normalize', true, valfcn_logical);
p.addParameter('mean_time', true, valfcn_logical);
p.addParameter('bin_size', 0, valfcn_posint);
p.addParameter('n_iter', 400, valfcn_posint);
p.addParameter('random_state', 0, valfcn_posint);
p.addParameter('kmeans_init_max', 20000, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('kmeans_max_samples', Inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('k', 0, valfcn_posint);
p.addParameter('cluster_method', 'hierarchical', @(x) any(strcmpi(x, {'hierarchical', 'contiguous'})));
p.addParameter('colors', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('clim', [0 1.5], @(x) validateattributes(x, {'numeric'}, {'numel', 2, 'increasing'}));

% Special key/value pairs that we have potentially set with optional keywords
p.addParameter('doplot', doplot, valfcn_logical);
p.addParameter('verbose', verbose, valfcn_logical);
p.addParameter('doverbose', verbose, valfcn_logical);

% process inputs
p.parse(X, varargin{:});

ARGS = p.Results;
ARGS.verbose = ARGS.verbose && ARGS.doverbose;
ARGS.cluster_method = lower(char(ARGS.cluster_method));
ARGS = rmfield(ARGS, {'X', 'doverbose'});

% END INPUT PARSER

if ARGS.n_clusters >= 200
    error('canlab_sort_rastermap:nClusters', 'n_clusters must be < 200 (the cluster sort is NP-hard).');
end

t_start = tic;
timing = struct();
vprint = @(varargin) verbose_print(ARGS.verbose, varargin{:});

% Reproducible seeding without disturbing the caller's RNG stream
rng_state = rng;
restore_rng = onCleanup(@() rng(rng_state)); %#ok<NASGU>
rng(ARGS.random_state, 'twister');

% -------------------------------------------------------------------------
% 1. Normalize
% -------------------------------------------------------------------------

[n_rows, n_cols] = size(X);
X = single(X);

if ARGS.time_bin > 1
    X = bin_columns(X, ARGS.time_bin);
end

igood = all(isfinite(X), 2);

if ARGS.normalize
    X = X - mean(X, 2);
    stdx = sqrt(mean(X.^2, 2));                % population SD, as numpy
    igood = igood & stdx > 0;
    X = X ./ stdx;
end

if ARGS.mean_time
    xm = mean(X(igood, :), 1);                 % across-row mean at each column
    xm = xm ./ norm(xm);
    if all(isfinite(xm))
        X = X - (X * xm') * xm;
    end
end

if ~ARGS.normalize
    igood = igood & std(X, 0, 2) > 0;
end

Xg = X(igood, :);
n_good = size(Xg, 1);

if n_good < 2
    error('canlab_sort_rastermap:tooFewRows', 'Fewer than 2 valid (finite, non-constant) rows to sort.');
end

timing.normalize = toc(t_start);
vprint('sorting %d valid rows x %d columns (%d rows excluded)\n', n_good, size(X, 2), n_rows - n_good);

% -------------------------------------------------------------------------
% 2. PCA: Usv = U*S (rows x PCs), Vsv = V*S (columns x PCs)
% -------------------------------------------------------------------------

t0 = tic;
n_PCs = min([ARGS.n_PCs, n_good, size(Xg, 2)]);
Usv = top_pcs(Xg, n_PCs);
sv = sqrt(sum(Usv.^2, 1));
keep_pc = sv > eps(single(max(sv))) * 10;
Usv = Usv(:, keep_pc); sv = sv(keep_pc);
n_PCs = size(Usv, 2);
Vsv = Xg' * (Usv ./ sv);                       % columns x PCs (= V*S)
timing.pca = toc(t0);
vprint('%d PCs computed, %.2f sec\n', n_PCs, timing.pca);

% -------------------------------------------------------------------------
% 3. Cluster rows (scaled k-means)
% -------------------------------------------------------------------------

t0 = tic;
Usv = double(Usv);

if n_good <= 50
    warning('canlab_sort_rastermap:fewRows', 'Data have <= 50 valid rows; skipping clustering and sorting rows directly.');
    U_nodes = Usv ./ max(sqrt(sum(Usv.^2, 2)), eps);
else
    n_clusters = min(floor(n_good / 2), ARGS.n_clusters);
    if isfinite(ARGS.kmeans_max_samples) && n_good > ARGS.kmeans_max_samples
        isub = randperm(n_good, ARGS.kmeans_max_samples);
        U_nodes = scaled_kmeans(Usv(isub, :), n_clusters, ARGS.kmeans_init_max);
    else
        U_nodes = scaled_kmeans(Usv, n_clusters, ARGS.kmeans_init_max);
    end
end
timing.kmeans = toc(t0);
vprint('%d clusters computed, %.2f sec\n', size(U_nodes, 1), timing.kmeans);

% Cluster similarity matrix
if ARGS.time_lag_window > 0
    cc = compute_cc_tdelay(double(Vsv) ./ sv, U_nodes, ARGS.time_lag_window, ARGS.symmetric);
else
    cc = U_nodes * U_nodes';
end

% -------------------------------------------------------------------------
% 4. Sort clusters (matrix matching / traveling salesman)
% -------------------------------------------------------------------------

t0 = tic;
[cc, inds] = tsp_sort(cc, ARGS.locality, ARGS.n_iter);
U_nodes = U_nodes(inds, :);
[~, embedding_clust] = max(Usv * U_nodes', [], 2);
timing.sort = toc(t0);
vprint('clusters sorted, %.2f sec\n', timing.sort);

% -------------------------------------------------------------------------
% 5. Upsample and embed rows
% -------------------------------------------------------------------------

t0 = tic;
nc = size(U_nodes, 1);
if ARGS.grid_upsample > 0 && nc >= 3        % local linear fits need >= 3 nodes
    gu = max(2, ARGS.grid_upsample);
    n_neighbors = max(min(8, nc - 1), floor(nc / 5));
    e_neighbor = max(1, min(1, n_neighbors - 1));
    [Y, gnode, Xrec, xy] = grid_upsampling(Usv, U_nodes, gu, n_neighbors, e_neighbor);
else
    % no upsampling: order rows by cluster, ties in original order
    [~, ord] = sort(embedding_clust);
    Y = zeros(n_good, 1); Y(ord) = 0:n_good - 1;
    gnode = embedding_clust;
    Xrec = U_nodes; xy = (0:nc - 1)';
end
timing.upsample = toc(t0);
vprint('clusters upsampled, %.2f sec\n', timing.upsample);

embedding = nan(n_rows, 1);
embedding(igood) = Y;
[~, isort] = sort(embedding);                  % NaNs (invalid rows) sort last

% -------------------------------------------------------------------------
% Optional: cluster rows (CANlab addition)
% -------------------------------------------------------------------------

cluster_labels = [];
linkage_tree = [];
if ARGS.k > 0
    t0 = tic;
    [lab_good, linkage_tree] = cluster_rows(Xrec, gnode, ARGS.k, ARGS.cluster_method);
    cluster_labels = nan(n_rows, 1);
    cluster_labels(igood) = lab_good;
    % number clusters in order of first appearance from the top of the display
    seen = cluster_labels(isort); seen = seen(~isnan(seen));
    [~, first] = unique(seen, 'stable');
    order = seen(sort(first));
    relab = zeros(max(order), 1); relab(order) = 1:numel(order);
    cluster_labels(igood) = relab(cluster_labels(igood));
    timing.cluster = toc(t0);
    vprint('%d row clusters (%s), %.2f sec\n', numel(order), ARGS.cluster_method, timing.cluster);
end

% -------------------------------------------------------------------------
% Display matrix: sorted, row-binned, z-scored
% -------------------------------------------------------------------------

bin_size = ARGS.bin_size;
if bin_size == 0 || n_good < bin_size || (bin_size == 50 && n_good < 1000)
    bin_size = max(1, floor(n_good / 500));
end
isort_good = isort(1:n_good);
X_embedding = bin_rows(X(isort_good, :), bin_size);
X_embedding = zscore_rows(X_embedding);

bin_labels = [];
if ~isempty(cluster_labels)
    lab_sorted = cluster_labels(isort_good);
    nb = floor(n_good / bin_size);
    bin_labels = mode(reshape(lab_sorted(1:nb * bin_size), bin_size, nb), 1)';
end

% -------------------------------------------------------------------------
% Collect outputs
% -------------------------------------------------------------------------

Vnorm = sqrt(sum(double(Vsv).^2, 2));
OUT = struct();
OUT.embedding = embedding;
OUT.isort = isort;
OUT.igood = igood;
OUT.cluster_labels = cluster_labels;
OUT.embedding_clust = embedding_clust;
OUT.X_embedding = X_embedding;
OUT.bin_size = bin_size;
OUT.bin_labels = bin_labels;
OUT.U_nodes = U_nodes;
OUT.X_nodes = U_nodes * (double(Vsv) ./ Vnorm)';
OUT.cc = cc;
OUT.U_upsampled = Xrec;
OUT.grid_positions = xy;
OUT.linkage_tree = linkage_tree;
OUT.settings = ARGS;
OUT.settings.n_PCs_used = n_PCs;
OUT.settings.n_clusters_used = nc;
timing.total = toc(t_start);
OUT.timing = timing;
OUT.fig_handle = [];

vprint('rastermap complete, %.2f sec\n', timing.total);

if ARGS.doplot
    OUT.fig_handle = plot_rastermap(OUT, ARGS);
end

end % main function


% =========================================================================
%
% Sub-functions
%
% =========================================================================

function verbose_print(doit, varargin)
if doit, fprintf(['canlab_sort_rastermap: ' varargin{1}], varargin{2:end}); end
end


function Xb = bin_columns(X, b)
% mean over bins of b adjacent columns (drops the remainder, as bin1d)
nb = floor(size(X, 2) / b);
Xb = squeeze(mean(reshape(X(:, 1:nb * b), size(X, 1), b, nb), 2));
if nb == 1, Xb = Xb(:); end
end


function Xb = bin_rows(X, b)
% mean over bins of b adjacent rows (drops the remainder, as bin1d)
if b <= 1, Xb = double(X); return, end
nb = floor(size(X, 1) / b);
Xb = squeeze(mean(reshape(double(X(1:nb * b, :)), b, nb, size(X, 2)), 1));
if nb == 1, Xb = Xb(:)'; end
end


function Z = zscore_rows(X)
X = X - mean(X, 2);
s = sqrt(mean(X.^2, 2)); s(s == 0) = 1;
Z = X ./ s;
end


function Usv = top_pcs(X, k)
% U*S for the top k principal axes of X (rows x cols), exact, via the
% eigendecomposition of the smaller Gram matrix.
[n, m] = size(X);
if m <= n
    C = double(X' * X);  C = (C + C') / 2;
    [V, D] = eig(C);
    [~, o] = sort(diag(D), 'descend');
    V = V(:, o(1:k));
    Usv = X * single(V);
else
    C = double(X * X');  C = (C + C') / 2;
    [U, D] = eig(C);
    [d, o] = sort(diag(D), 'descend');
    U = U(:, o(1:k));
    Usv = single(U .* sqrt(max(d(1:k), 0))');
end
% sign convention: largest-magnitude loading positive (deterministic output)
[~, imax] = max(abs(Usv), [], 1);
sgn = sign(Usv(sub2ind(size(Usv), imax, 1:size(Usv, 2))));
sgn(sgn == 0) = 1;
Usv = Usv .* sgn;
end


function X_nodes = scaled_kmeans(X, n_clusters, init_max)
% k-means with correlation distance and k-means++ seeding (cluster.py)

n_iter = 50;
n_samples = size(X, 1);

if n_samples > init_max
    Xi = X(randperm(n_samples, init_max), :);
else
    Xi = X;
end
X_nodes = scaled_kmeans_init(Xi, n_clusters, 100);
X_nodes = X_nodes ./ sqrt(1e-4 + sum(X_nodes.^2, 2));

Xs = single(X);
for j = 1:n_iter
    cc = Xs * single(X_nodes');
    [mx, imax] = max(cc, [], 2);
    % each row contributes to its best cluster, weighted by its similarity
    W = sparse(1:n_samples, double(imax), double(mx), n_samples, n_clusters);
    X_nodes = full(W' * X);
    X_nodes = X_nodes ./ sqrt(1e-10 + sum(X_nodes.^2, 2));
end

X_nodes = X_nodes(sqrt(sum(X_nodes.^2, 2)) > 0, :);
[~, o] = sort(X_nodes(:, 1));
X_nodes = X_nodes(o, :);

if size(X_nodes, 1) < floor(n_clusters / 2)
    warning('canlab_sort_rastermap:kmeans', 'Found fewer than half the requested clusters; consider reducing n_clusters.');
end
end


function X_nodes = scaled_kmeans_init(X, n_clusters, n_local_trials)
% k-means++ seeding with correlation distance (_scaled_kmeans_init)
n_samples = size(X, 1);
n_local_trials = min(n_samples - 1, n_local_trials);
X_nodes = zeros(n_clusters, size(X, 2));

X_norm = sqrt(sum(X.^2, 2));
cid = randi(n_samples);
X_nodes(1, :) = X(cid, :);
closest = 1 - (X * X_nodes(1, :)') ./ (X_norm * norm(X_nodes(1, :)));
pot = sum(closest);

for c = 2:n_clusters
    rand_vals = rand(n_local_trials, 1) * pot;
    cdf = cumsum(closest);
    cand = min(arrayfun(@(v) find(cdf >= v, 1, 'first'), rand_vals, 'UniformOutput', true), n_samples);
    Xc = X(cand, :);
    D = 1 - (X * Xc') ./ (X_norm * sqrt(sum(Xc.^2, 2))');
    newd = min(closest, D);
    pots = sum(newd, 1);
    [pot, best] = min(pots);
    X_nodes(c, :) = Xc(best, :);
    closest = newd(:, best);
end
end


function cc = compute_cc_tdelay(V, U_nodes, w, symmetric)
% max cross-correlation of cluster time courses over time lags (cluster.py)
Xn = U_nodes * V';
Xn = zscore_rows(Xn);
nt = size(Xn, 2);
cc = -Inf(size(Xn, 1));
for tshift = -w * symmetric : w
    if tshift < 0
        c = (Xn(:, 1:nt + tshift) * Xn(:, 1 - tshift:end)') / (nt - tshift);
    else
        c = (Xn(:, 1 + tshift:end) * Xn(:, 1:nt - tshift)') / (nt - tshift);
    end
    cc = max(cc, c);
end
end


% -------------------------------------------------------------------------
% Cluster sort: maximize sum(cc(p,p) .* BBt) over orderings p (sort.py)
% -------------------------------------------------------------------------

function [cc, inds] = tsp_sort(cc, locality, n_iter)
n = size(cc, 1);
inds = (1:n)';
if n < 3, return, end

b = bbt_profile(n, locality);                  % b(d+1) = template at distance d
n_skip = max(1, floor(n / 30));

[cc, inds2] = tsp_pass(cc, b, n_skip, n_iter);
inds = inds(inds2);
if n_skip > 1
    [cc, inds2] = tsp_pass(cc, b, 1, n_iter);
    inds = inds(inds2);
end
end


function b = bbt_profile(n, locality)
% Template value as a function of distance d = 0..n between positions in the
% ordering (compute_BBt with xi = (0:n-1)/n), including the normalization by
% the mean over the full n x n matrix. Diagonal (d = 0) is 0.
epsl = 1e-3;
d = (0:n)';
bl = -log(epsl + d / n); bl(1) = 0;
% mean over the full n x n matrix: distance d occurs (n-d)*2 times (d>0)
cnt = [n; 2 * (n - (1:n - 1)')];
bl = bl / (sum(bl(1:n) .* cnt) / n^2);
if locality > 0
    bm = exp(-(d / n).^2 / (2 * (0.5 / n)^2)); bm(1) = 0;
    bm = bm / (sum(bm(1:n) .* cnt) / n^2);
    b = bl * (1 - locality) + bm * locality;
else
    b = bl;
end
end


function [G, inds] = tsp_pass(G, b, n_skip, n_iter)
% One pass of tsp_fast: greedy best segment move per iteration, testing
% segment lengths in pairs (1-2, 3-4, ...) and stopping at the first pair
% that yields an improvement > 1e-3.

n = size(G, 1);
inds = (1:n)';
n_test = floor(n / n_skip);
B = triu(toeplitz(b(1:n)));                     % B(i,j) = b(j-i), j >= i
corr_orig = sum(G .* B, 'all');
n_pairs = floor(n / 2);
cache = cell(n, 1);

for k = 0:n_iter - 1
    improved = false;
    for tl = 1:n_pairs
        lens = [2 * tl - 1, 2 * tl];
        change = -Inf(n_test, n, 2);
        isrev = false(n_test, n, 2);
        for li = 1:2
            L = lens(li);
            if L > n, continue, end
            m = n - L;
            q = mod(k, n_skip) + (0:n_test - 1) * n_skip;   % candidate new start positions
            qok = q <= m;
            if isempty(cache{L}), cache{L} = move_cache(n, L, b); end
            f = score_moves(G, L, b, cache{L}, false);       % (m+1) x (m+1): q x i0
            r = false(size(f));
            if L > 1
                fr = score_moves(G, L, b, cache{L}, true);
                r = fr > f;
                f(r) = fr(r);
            end
            c = -Inf(n_test, n);
            c(qok, 1:m + 1) = f(q(qok) + 1, :) - corr_orig;
            rr = false(n_test, n);
            rr(qok, 1:m + 1) = r(q(qok) + 1, :);
            change(:, :, li) = c;
            isrev(:, :, li) = rr;
        end
        [best, ix] = max(change(:));
        if best > 1e-3
            improved = true;
            [qi, i0, li] = ind2sub(size(change), ix);
            L = lens(li);
            q = mod(k, n_skip) + (qi - 1) * n_skip;
            i0 = i0 - 1;
            s = i0 + (1:L)';
            if isrev(ix), s = flipud(s); end
            rest = [1:i0, i0 + L + 1:n]';
            perm = [rest(1:q); s; rest(q + 1:end)];
            G = G(perm, perm);
            inds = inds(perm);
            corr_orig = sum(G .* B, 'all');
            break
        end
    end
    if ~improved, break, end
end
end


function C = move_cache(n, L, b)
% Index/weight arrays that depend only on n and the segment length L
m = n - L;
C = struct();
a = (1:L)';
% Part 1 (rest elements before the segment): H1(k,c) = sum_a Grs(k,a) b(c-1+a)
C.Bh1 = reshape(b(a + (1:m + 1)), L, m + 1); % b index (c-1+a)+1
% Part 2 (rest elements after the segment): H2(k,j) = sum_a Gsr(k,a) b(j+L-a)
C.Bh2 = reshape(b((1:m) + L - a + 1), L, m);
% anti-diagonal sums of H1 (k + c - 1 = q) and diagonal sums of H2 (k - j = q)
[kk, cc] = ndgrid(1:m, 1:m + 1);
qq = kk + cc - 1;
ok = qq <= m;
C.A1 = sparse(qq(ok) + 1, find(ok), 1, m + 1, m * (m + 1));
[kk, jj] = ndgrid(1:m, 1:m);
qq = kk - jj;
ok = qq >= 0;
C.A2 = sparse(qq(ok) + 1, find(ok), 1, m + 1, m * m);
% remainder pairs: template at distance l-k, and change if the segment is
% inserted between them (distance grows by L)
D = (1:m)' - (1:m);                             % k - l
D = -D;                                         % l - k
up = D > 0;
Bm = zeros(m); Bm(up) = b(D(up) + 1);
Em = zeros(m); Em(up) = b(D(up) + L + 1) - b(D(up) + 1);
C.Bm = Bm; C.Em = Em;
Ds = (1:L) - (1:L)';                            % b2 - a
Bs = zeros(L); Bs(Ds > 0) = b(Ds(Ds > 0) + 1);
C.Bs = Bs;
% linear indices of T(q, q+1), q = 1..m-1, for the remainder-split term
q = (1:m - 1)';
C.Tidx = q + q * m;
end


function f = score_moves(G, L, b, C, reversed) %#ok<INUSL>
% Score sum(Gnew .* B) for moving the segment at positions i0+1..i0+L to new
% start position q, for all i0 = 0..m and q = 0..m (m = n-L). Returns
% f(q+1, i0+1). If reversed, the segment is also reversed.
n = size(G, 1);
m = n - L;
P = m + 1;
i0 = 0:m;
Sidx = (1:L)' + i0;                             % L x P
if reversed, Sidx = flipud(Sidx); end
Ridx = (1:m)' + L * ((1:m)' > i0);              % m x P

Rk = reshape(Ridx, m, 1, P);
Sa = reshape(Sidx, 1, L, P);
Grs = G(Rk + (Sa - 1) * n);                     % G(r_k, s_a): rest before segment
Gsr = G(Sa + (Rk - 1) * n);                     % G(s_a, r_k): segment before rest

H1 = pagemtimes(Grs, C.Bh1);                    % m x (m+1) x P
H2 = pagemtimes(Gsr, C.Bh2);                    % m x m x P
part1 = C.A1 * reshape(H1, [], P);              % (m+1) x P
part2 = C.A2 * reshape(H2, [], P);

GR = G(Rk + (reshape(Ridx, 1, m, P) - 1) * n);  % m x m x P (rest x rest)
base = reshape(sum(GR .* C.Bm, [1 2]), 1, P);
E = GR .* C.Em;
T = cumsum(cumsum(E, 1), 2, 'reverse');         % T(q,l) = sum_{k<=q, l'>=l} E
split = zeros(m + 1, P);
if m > 1
    split(2:m, :) = T(C.Tidx + (0:P - 1) * m * m);
end

GS = G(reshape(Sidx, L, 1, P) + (Sa - 1) * n);  % L x L x P (segment, oriented)
fss = reshape(sum(GS .* C.Bs, [1 2]), 1, P);

f = base + fss + split + part1 + part2;
end


% -------------------------------------------------------------------------
% Upsampling (upsample.py)
% -------------------------------------------------------------------------

function [Y, imax, Xrec, xy] = grid_upsampling(X, X_nodes, gu, n_neighbors, e_neighbor)
nc = size(X_nodes, 1);
Y_nodes = (0:nc - 1)';
xy = (0:1 / gu:nc - 1 / gu)';                   % grid positions
xy = xy(1:nc * gu);
nX = numel(xy);
ds = (xy - Y_nodes').^2;                        % nX x nc
[dss, isrt] = sort(ds, 2);
Xrec = zeros(nX, size(X_nodes, 2));
for j = 1:nX
    ineigh = isrt(j, 1:n_neighbors);
    dists = dss(j, 1:n_neighbors);
    w = exp(-dists / dists(e_neighbor + 1))';
    M = X_nodes(ineigh, :);
    N = [Y_nodes(ineigh), ones(n_neighbors, 1)];
    A = N' * (N .* w);                          % weighted local linear fit:
    Ainv = [A(2, 2), -A(1, 2); -A(2, 1), A(1, 1)] / (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1));
    R = Ainv * (N' * (M .* w));                 % (explicit 2x2 inverse: mldivide with
                                                %  many right-hand sides is ~1000x slower here)
    Xrec(j, :) = xy(j) * R(1, :) + R(2, :);
end
Xrec = Xrec ./ sqrt(sum(Xrec.^2, 2));

% assign each row to the grid node it correlates with best (in chunks)
Xz = zscore_rows(double(X));
n = size(Xz, 1);
imax = zeros(n, 1);
chunk = max(1, floor(2e7 / nX));
Xr = single(Xrec);
for i0 = 1:chunk:n
    ii = i0:min(n, i0 + chunk - 1);
    c = max(0, Xr * single(Xz(ii, :)'));
    [~, imax(ii)] = max(c, [], 1);
end
Y = xy(imax);
end


% -------------------------------------------------------------------------
% Row clustering (CANlab addition)
% -------------------------------------------------------------------------

function [lab, Z] = cluster_rows(Xrec, gnode, k, method)
% Cluster the occupied grid nodes and give each row its node's label.
occ = unique(gnode);
k = min(k, numel(occ));
Z = [];
switch method
    case 'hierarchical'
        if numel(occ) > 1
            Z = linkage(Xrec(occ, :), 'ward');
            lab_occ = cluster(Z, 'maxclust', k);
        else
            lab_occ = 1;
        end
    case 'contiguous'
        counts = accumarray(gnode, 1, [size(Xrec, 1), 1]);
        lab_occ = contiguous_ward(Xrec(occ, :), counts(occ), k);
end
map = zeros(size(Xrec, 1), 1);
map(occ) = lab_occ;
lab = map(gnode);
end


function lab = contiguous_ward(V, w, k)
% Agglomerate adjacent segments (in row order of V) with the Ward criterion,
% weighting each node by w, until k segments remain.
n = size(V, 1);
starts = (1:n)';
sz = w(:);
S = V .* sz;                                    % weighted sums
nseg = n;
while nseg > k
    mu = S ./ sz;
    dmu = sum((mu(1:end - 1, :) - mu(2:end, :)).^2, 2);
    cost = (sz(1:end - 1) .* sz(2:end)) ./ (sz(1:end - 1) + sz(2:end)) .* dmu;
    [~, i] = min(cost);
    S(i, :) = S(i, :) + S(i + 1, :);  S(i + 1, :) = [];
    sz(i) = sz(i) + sz(i + 1);        sz(i + 1) = [];
    starts(i + 1) = [];
    nseg = nseg - 1;
end
lab = zeros(n, 1);
lab(starts) = 1;
lab = cumsum(lab);
end


% -------------------------------------------------------------------------
% Plot (CANlab addition)
% -------------------------------------------------------------------------

function fh = plot_rastermap(OUT, ARGS)
Xe = OUT.X_embedding;
clim = ARGS.clim;
fh = figure('Color', 'w', 'Name', 'canlab_sort_rastermap');

if OUT.bin_size > 1
    rowlabel = sprintf('Sorted rows (bins of %d)', OUT.bin_size);
else
    rowlabel = 'Sorted rows';
end

if isempty(OUT.bin_labels)
    ax = axes(fh);
    imagesc(ax, Xe, clim);
    colormap(ax, flipud(gray(256)));
    ylabel(ax, rowlabel);
else
    k = max(OUT.bin_labels);
    cols = parse_colors(ARGS.colors, k);
    % colour each bin's row in its cluster colour: white (low) -> colour (high)
    v = min(max((Xe - clim(1)) / diff(clim), 0), 1);
    C = cols(OUT.bin_labels, :);                    % n_bins x 3
    rgb = 1 - v .* reshape(1 - C, [], 1, 3);
    ax = axes(fh, 'Position', [0.15 0.11 0.80 0.815]);
    image(ax, rgb);
    set(ax, 'YTick', []);
    % cluster colour strip along the rows, labelled with cluster numbers
    axs = axes(fh, 'Position', [0.10 0.11 0.035 0.815]);
    image(axs, reshape(C, [], 1, 3));
    set(axs, 'XTick', [], 'YDir', 'reverse', 'TickLength', [0 0]);
    [~, firstbin] = unique(OUT.bin_labels, 'stable');
    lastbin = [firstbin(2:end) - 1; numel(OUT.bin_labels)];
    set(axs, 'YTick', (firstbin + lastbin) / 2, 'YTickLabel', OUT.bin_labels(firstbin));
    axis(axs, 'tight');
    ylabel(axs, rowlabel);
    linkaxes([ax axs], 'y');
end

set(ax, 'YDir', 'reverse');
axis(ax, 'tight');
xlabel(ax, 'Columns (e.g., time / images)');
title(ax, 'Rastermap-sorted data');
try canlab_hide_axes_toolbar(fh); catch, end
end


function cols = parse_colors(colors, k)
% colors: [] (seaborn), [n x 3] matrix, colormap name, or cell of [r g b]
if isempty(colors)
    [~, cols] = seaborn_colors(k);
elseif iscell(colors)
    cols = cell2mat(cellfun(@(c) c(:)', colors(:), 'UniformOutput', false));
elseif ischar(colors) || isstring(colors)
    name = char(colors);
    if strcmpi(name, 'seaborn')
        [~, cols] = seaborn_colors(k);
    else
        cols = feval(name, k);
    end
else
    cols = colors;
end
if size(cols, 2) ~= 3
    error('canlab_sort_rastermap:colors', 'Colors must be RGB triplets (n x 3).');
end
if size(cols, 1) < k                             % recycle colours if too few
    cols = cols(mod(0:k - 1, size(cols, 1)) + 1, :);
end
cols = cols(1:k, :);
end
