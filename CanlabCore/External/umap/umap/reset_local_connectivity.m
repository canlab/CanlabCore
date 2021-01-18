function simplicial_set = reset_local_connectivity(simplicial_set)
%RESET_LOCAL_CONNECTIVITY Reset the local connectivity requirement -- each
% data sample should have complete confidence in at least one 1-simplex in
% the simplicial set. We can enforce this by locally rescaling confidences,
% and then remerging the different local simplicial sets together.
% 
% simplicial_set = RESET_LOCAL_CONNECTIVITY(simplicial_set)
%
% Parameters
% ----------
% simplicial_set: sparse matrix
%     The simplicial set for which to recalculate with respect to local
%     connectivity.
% 
% Returns
% -------
% simplicial_set: sparse matrix
%     The recalculated simplicial set, now with the local connectivity
%     assumption restored.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
    
    n_cols = size(simplicial_set, 2);

    divisor = max(simplicial_set,[],2);
    if ~issparse(simplicial_set)
        simplicial_set = simplicial_set./repmat(divisor, [1 n_cols]);
    else
        simplicial_set = divide_sparse(simplicial_set, divisor);
    end

    transpose = simplicial_set';
    prod_matrix = simplicial_set.*transpose;
    simplicial_set = simplicial_set + transpose - prod_matrix;