function result = general_simplicial_set_intersection(simplicial_set1, simplicial_set2, weight)
%GENERAL_SIMPLICIAL_SET_INTERSECTION Combine two fuzzy simplicial sets via
% intersection.
%
% result = GENERAL_SIMPLICIAL_SET_INTERSECTION(simplicial_set1, simplicial_set2, weight)
% 
% Parameters
% ----------
% simplicial_set1: sparse matrix
%     The first input fuzzy simplicial set.
% 
% simplicial_set2: sparse matrix
%     The second input fuzzy simplicial set.
% 
% weight: double (optional, default 0.5)
%     The weight assigned to simplicial_set1 when performing the
%     intersection.
% 
% Returns
% -------
% result: sparse matrix
%     The resulting intersected fuzzy simplicial set.
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

    if nargin < 3
        weight = 0.5;
    end

    result = sparse(simplicial_set1 + simplicial_set2);
    left = sparse(simplicial_set1);
    right = sparse(simplicial_set2);
    
    [n_row, n_col] = size(result);
    
    [result_row, result_col, result_data] = find(result);

    %Not clear if I have to do this for MATLAB
    result_values = general_sset_intersection(left, right, result_row, result_col,...
        result_data, weight);
    
    result = sparse(result_row, result_col, result_values, n_row, n_col);
    
end