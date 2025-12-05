function X = intercept(X, meth)
% Intercept-related functions for working with design matrices, etc.
%
% Return which columns are intercept(s), if any
% ::
%
%    wh = intercept(X, 'which')
%
% Remove an intercept, if there is one
% ::
%
%    X = intercept(X, 'remove');
%
% Add an intercept to the end
% ::
%
%    X = intercept(X, 'add');
%
% Ensure that the intercept is at the end, moving or adding it as necessary
% ::
%
%    X = intercept(X, 'end');
%
%   General usability edits and edge-case fixes by Michael Sun, Ph.D. on
%   12/05/2025

if nargin < 2
    meth = 'which';
end

tol = 1e-8;  % tolerance for floating point

% --- Find intercept columns: columns that are (approximately) all ones ---
[nObs, nCols] = size(X);

if nCols == 0 || nObs == 0
    % Empty design matrix or no observations
    wh_is_intercept = [];
else
    % A column is an intercept if all entries are ~1
    is_intercept_col = all(abs(X - 1) < tol, 1);
    wh_is_intercept  = find(is_intercept_col);
end

switch lower(meth)
    
    case 'which'
        X = wh_is_intercept;
        
    case 'remove'
        if ~isempty(wh_is_intercept)
            X(:, wh_is_intercept) = [];
        end
        
    case 'add'
        X(:, end+1) = 1;
        
    case 'end'
        X = intercept(X, 'remove');
        X = intercept(X, 'add');
        
    otherwise
        error('Unknown method');
        
end

end
        
