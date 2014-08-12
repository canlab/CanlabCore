% function [indic, xlevels] = condf2indic(X)
%
% Create N x k indicator matrix of ones/zeros from N x 1 list of numeric values (X)
%
% indic is returned as single precision type, because it can then be used
% as a design matrix in a GLM
% xlevels are the values of X corresponding to columns of indic
%
% tor wager, nov 2007

function [indic, xlevels] = condf2indic(X)

    [xlevels, ind1, indx] = unique(X);

    indic = single(false(length(indx), length(xlevels)));

    for i = 1:length(xlevels)
        indic(indx == i, i) = 1;
    end

end