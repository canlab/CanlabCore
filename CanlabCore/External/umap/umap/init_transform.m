function new_embedding = init_transform(indices, weights, embedding)
%INIT_TRANSFORM Given indices and weights and an original embedding
% initialize the positions of new points relative to the
% indices and weights (of their neighbors in the source data).
%
% new_embedding = INIT_TRANSFORM(indices, weights, embedding)
%
% Parameters
% ----------
% indices: array of size (n_new_samples, n_neighbors)
%     The indices of the neighbors of each new sample
% 
% weights: array of size (n_new_samples, n_neighbors)
%     The membership strengths of associated 1-simplices
%     for each of the new samples.
% 
% embedding: array of size (n_samples, dim)
%     The original embedding of the source data.
% 
% Returns
% -------
% new_embedding: array of size (n_new_samples, dim)
%     An initial embedding of the new sample points.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

    
    [A1,A2,A3] = ndgrid(1:size(indices,1), 1:size(embedding,2), 1:size(indices,2));
    array = arrayfun(@(x,y,z) weights(x,z)*embedding(indices(x,z), y), A1, A2, A3);

    new_embedding = sum(array, 3);