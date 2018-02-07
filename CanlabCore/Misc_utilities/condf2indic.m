function [indic, xlevels] = condf2indic(X, varargin)
% Create N x k indicator matrix of ones/zeros from N x 1 list of numeric
% values (X)
%
% :Usage:
% ::
%
%     [indic, xlevels] = condf2indic(X, ['integers'])
%
% :Inputs:
%
%   **'integers':**
%       Treat X as integer vector. Check that all are integers, and return
%       one column per integer in order, excluding zeros. Insert NaNs in
%       columns whose corresponding integers are missing.  Returns an n x k
%       matrix indic, where n is size(X, 1) and k is max(X).
%
% :Outputs:
%
%   **indic:**
%        is returned as single precision type, because it can then be
%        used as a design matrix in a GLM
%
%   **xlevels:**
%        are the values of X corresponding to columns of indic
%
% ..
%    tor wager, nov 2007. Edited Feb 2018 to expand functionality to
%    'integers' case, for use with atlas object.
% ..

X = double(X); % for objects/special types

[xlevels, ind1, indx] = unique(X);

indic = single(false(length(indx), length(xlevels)));

for i = 1:length(xlevels)
    indic(indx == i, i) = 1;
end

if any(strcmp(varargin, 'integers'))
    
    % Check integers
    u = unique(X);
    if ~all(u == round(u))
        warning('Some condition function values are not integers.');
    end
    
    % Remove zero-valued integer
    wh = find(xlevels == 0);
    xlevels(wh) = [];
    indic(:, wh) = [];
    
    % if codes are 1...n, xlevels is index values...but if there are
    % missing integers, need to insert into indic
    
    is_missing = true(max(X), 1);
    is_missing(xlevels) = false;
    
    % Insert NaNs for any missing index values.
    
    indic = naninsert(is_missing, indic')';
    
end % integers

end % main function
