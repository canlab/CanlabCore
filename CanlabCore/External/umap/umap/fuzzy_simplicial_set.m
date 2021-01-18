function fuzzy_simplicial_set = fuzzy_simplicial_set(X, n_neighbors, metric, varargin)
%FUZZY_SIMPLICIAL_SET Given a set of data X, a neighborhood size, and a
% measure of distance compute the fuzzy simplicial set (here represented as
% a fuzzy graph in the form of a sparse matrix) associated to the data.
% This is done by locally approximating geodesic distance at each point,
% creating a fuzzy simplicial set for each such point, and then combining
% all the local fuzzy simplicial sets into a global one via a fuzzy union.
%
% fuzzy_simplicial_set = FUZZY_SIMPLICIAL_SET(X, n_neighbors, metric)
% 
% Parameters
% ----------
% X: array of size (n_samples, n_features)
%     The data to be modelled as a fuzzy simplicial set.
% 
% n_neighbors: double
%     The number of neighbors to use to approximate geodesic distance.
%     Larger numbers induce more global estimates of the manifold that can
%     miss finer detail, while smaller values will focus on fine manifold
%     structure to the detriment of the larger picture.
% 
% metric: string or function (optional, default 'euclidean')
%     The metric to use to compute distances in high dimensional space.
%     If a string is passed it must match a valid predefined metric. For now,
%     valid string metrics include:
%         * euclidean (or l2)
%         * manhattan (or l1)
%         * chebyshev (or linf)
%         * correlation
%         * cosine
%         * hamming
%         * jaccard
%         * mahalanobis
%         * minkowski
%         * seuclidean
% 
% 
% knn_indices: array of size (n_samples, n_neighbors) (optional)
%     If the k-nearest neighbors of each point has already been calculated
%     you can pass them in here to save computation time. This should be
%     an array with the indices of the k-nearest neighbors as a row for
%     each data point.
% 
% knn_dists: array of size (n_samples, n_neighbors) (optional)
%     If the k-nearest neighbors of each point has already been calculated
%     you can pass them in here to save computation time. This should be
%     an array with the distances of the k-nearest neighbors as a row for
%     each data point.
% 
% set_op_mix_ratio: double (optional, default 1)
%     Interpolate between (fuzzy) union and intersection as the set operation
%     used to combine local fuzzy simplicial sets to obtain a global fuzzy
%     simplicial sets. Both fuzzy set operations use the product t-norm.
%     The value of this parameter should be between 0 and 1; a value of
%     1 will use a pure fuzzy union, while 0 will use a pure fuzzy
%     intersection.
% 
% local_connectivity: double (optional, default 1)
%     The local connectivity required -- i.e. the number of nearest
%     neighbors that should be assumed to be connected at a local level.
%     The higher this value the more connected the manifold becomes
%     locally. In practice this should be not more than the local intrinsic
%     dimension of the manifold.
% 
% umap is the umap instance
%
% Returns
% -------
% fuzzy_simplicial_set: sparse matrix of size (n_samples, n_samples)
%     A fuzzy simplicial set represented as a sparse matrix. The (i,
%     j) entry of the matrix represents the membership strength of the
%     1-simplex between the ith and jth sample points.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

    p=parseArguments();
    parse(p,varargin{:});
    args=p.Results; 
    local_connectivity = args.local_connectivity;
    set_op_mix_ratio = args.set_op_mix_ratio;
    knn_dists = args.knn_dists;
    knn_indices = args.knn_indices;
    
    X_rows = size(X, 1);
    
    if isempty(knn_indices) || isempty(knn_dists)
        if isempty(args.umap)
            [knn_indices, knn_dists] = nearest_neighbors(X, ...
                struct('K', n_neighbors, 'Distance', metric));
        else
            [knn_indices, knn_dists] = nearest_neighbors(X, args.umap);
        end
    end

    [sigmas, rhos] = smooth_knn_dist(knn_dists, n_neighbors, 'local_connectivity', local_connectivity);

    [rows, cols, vals] = compute_membership_strengths(knn_indices, knn_dists, sigmas, rhos);

    result = sparse(rows, cols, vals, X_rows, X_rows);

    transpose = result';

    prod_matrix = result .* transpose;

    fuzzy_simplicial_set = (set_op_mix_ratio * (result + transpose - prod_matrix) +...
        (1 - set_op_mix_ratio) * prod_matrix);
    
    function p=parseArguments(varargin)
        p = inputParser;
        addParameter(p,'local_connectivity', 1, @(x) isnumeric(x));
        addParameter(p,'set_op_mix_ratio', 1, @(x) isnumeric(x));
        addParameter(p,'knn_dists', []);
        addParameter(p,'knn_indices', []);
        addParameter(p,'umap', [], @(x) isnumeric(x));
    end
    
end
        