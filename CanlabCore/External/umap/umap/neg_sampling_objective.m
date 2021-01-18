function result = neg_sampling_objective(head_embedding, tail_embedding, head, tail, weights, a, b, negative_sample_rate, same_embedding)
%NEG_SAMPLING_OBJECTIVE Given a low-dimensional embedding and weights of
% 1-simplices of a high-dimensional simplicial complex, compute the
% associated negative sampling objective. This is the quantity that the
% original Python UMAP implementation is attempting to minimize by using
% the negative sampling technique during stochastic gradient descent (SGD).
% Note that this value is distinct from cross entropy as presented in the
% original UMAP paper.
% This calculation uses the modified smooth formula Phi for
% low-dimensional weight that is used in SGD.
% Because the computation time is O(n^2), we recommend only using this when
% n1*n2 <= 1e8.
%
% result = neg_sampling_objective(head_embedding, tail_embedding, head, tail, weights, a, b, negative_sample_rate, same_embedding)
%
% Parameters
% ----------
% n_neighbors: double (optional, default 15)
%     The size of local neighborhood (in terms of number of neighboring
%     sample points) used for manifold approximation. Larger values result
%     in more global views of the manifold, while smaller values result in
%     more local data being preserved. In general values should be in the
%     range 2 to 100.
% 
% n_components: integer (optional, default 2)
%     The dimension of the space to embed into. This defaults to 2 to
%     provide easy visualization, but can reasonably be set to any integer
%     value in the range 2 to 100.
% 
% metric: string or function (optional, default 'euclidean')
%     The metric to use to compute distances in high dimensional space. If
%     a string is passed, it must match a valid predefined metric. For now,
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
% n_epochs: integer (optional)
%     The number of training epochs to be used in optimizing the low
%     dimensional embedding. Larger values result in more accurate
%     embeddings. If 0, a value will be selected based on the size of the
%     input dataset (200 for large datasets, 500 for small).
% 
% learning_rate: double (optional, default 1)
%     The initial learning rate for the embedding optimization.
% 
% init: string (optional, default 'spectral')
%     How to initialize the low dimensional embedding. Options are:
%         * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
%         * 'random': assign initial embedding positions at random.
%         * An array of initial embedding positions.
% 
% min_dist: double (optional, default 0.1)
%     The effective minimum distance between embedded points. Smaller
%     values will result in a more clustered/clumped embedding where nearby
%     points on the manifold are drawn closer together, while larger values
%     will result on a more even dispersal of points. The value should be
%     set relative to the "spread" value, which determines the scale at
%     which embedded points will be spread out.
% 
% spread: double (optional, default 1)
%     The effective scale of embedded points. In combination with
%     "min_dist" this determines how clustered/clumped the embedded points
%     are.
% 
% set_op_mix_ratio: double (optional, default 1)
%     Interpolate between (fuzzy) union and intersection as the set
%     operation used to combine local fuzzy simplicial sets to obtain a
%     global fuzzy simplicial sets. Both fuzzy set operations use the
%     product t-norm. The value of this parameter should be between 0 and
%     1; a value of 1 will use a pure fuzzy union, while 0 will use a pure
%     fuzzy intersection.
% 
% local_connectivity: integer (optional, default 1)
%     The local connectivity required -- i.e. the number of nearest
%     neighbors that should be assumed to be connected at a local level.
%     The higher this value the more connected the manifold becomes
%     locally. In practice this should be not more than the local intrinsic
%     dimension of the manifold.
% 
% repulsion_strength: double (optional, default 1)
%     Weighting applied to negative samples in low dimensional embedding
%     optimization. Values higher than one will result in greater weight
%     being given to negative samples.
% 
% negative_sample_rate: integer (optional, default 5)
%     The number of negative samples to select per positive sample in the
%     optimization process. Increasing this value will result in greater
%     repulsive force being applied, greater optimization cost, but
%     slightly more accuracy.
% 
% transform_queue_size: double (optional, default 4)
%     For transform operations (embedding new points using a trained model)
%     this will control how aggressively to search for nearest neighbors.
%     Larger values will result in slower performance but more accurate
%     nearest neighbor evaluation.
% 
% a: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% b: double
%     Parameter of differentiable approximation of right adjoint functor.ad
% 
% Returns
% -------
% result: double
%     The value of the negative sampling objective.
%
% See also: CROSS_ENTROPY
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
    n1 = size(head_embedding, 1);
    n2 = size(tail_embedding, 1);
    if n1*n2 > 1e8
        error('HALTED: MATLAB usually freezes for embeddings this large.');
    end
    
    n_other_points = n2;
    if same_embedding
        n_other_points = n_other_points - 1;
    end

    full_dists = pdist2(head_embedding, tail_embedding);

    Phi = ones(size(full_dists))./(1 + a*(full_dists.^(2*b)));
    if same_embedding
        Phi = Phi - diag(diag(Phi));
    end
    gap_part = sum(log(1 - Phi), 2);

    Phi_summands = weights.*(log(Phi(sub2ind(size(Phi), head, tail))) + negative_sample_rate/n_other_points*gap_part(head));

    result = -sum(Phi_summands);
end