function result = cross_entropy(head_embedding, tail_embedding, head, tail, weights, a, b, same_embedding)
%CROSS_ENTROPY Given a distance for each 1-simplex in low-dimensional space
% and the original weights of the 1-simplices in high-dimensional space,
% compute the approximation to the cross-entropy between the two simplicial
% complexes. This calculation uses the modified smooth formula Phi for
% low-dimensional weight that is used in the stochastic gradient descent.
%
% result = cross_entropy(head_embedding, tail_embedding, head, tail, weights, a, b, same_embedding)
%
% Parameters
% ----------
% dists: array of size (n_1_simplices, 1)
%     The current distance between the two endpoints of the 1-simplex in
%     low-dimensional Euclidean space.
%
% weights: array of size (n_1_simplices, 1)
%     The original weights assigned to the 1-simplices in high-dimensional
%     space.
%
% a: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% b: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% Returns
% -------
% result: double
%     The total approximated cross entropy between the two simplicial complexes.
%
% See also: NEG_SAMPLING_OBJECTIVE
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

    full_dists = pdist2(head_embedding, tail_embedding);
    full_weights = full(sparse(head, tail, weights, n1, n2));
    if same_embedding
        full_weights = full_weights + eye(n1);
    end
    Phi = ones(size(full_weights))./(1 + a*(full_dists.^(2*b)));
    fw0 = full_weights == 0;
    fw1 = full_weights == 1;
    other = ~fw0 & ~fw1;
    Phi_summands = zeros(size(full_weights));
    Phi_summands(fw0) = log(1-Phi(fw0));
    Phi_summands(fw1) = log(Phi(fw1));
    Phi_summands(other) = full_weights(other).*log(Phi(other)) + (1-full_weights(other)).*log(1-Phi(other));

    result = -sum(sum(Phi_summands));
end