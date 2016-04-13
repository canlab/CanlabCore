function subsets = subset_indicator_matrix(n)
% Create a matrix whose rows contain indicators (1/0 values) for all
% possible subsets of n variables
%
% :Usage:
% ::
%
%     subsets = subset_indicator_matrix(n)
%
% Use this, for example, to create a matrix that tells you all possible
% combinations of regressors to include in a regression.
%
%
%    Tor Wager, March 07
%


    vals = cell(1,n);
    [vals{:}] = deal([1 0]);
    subsets = combvec(vals{:})';

end
