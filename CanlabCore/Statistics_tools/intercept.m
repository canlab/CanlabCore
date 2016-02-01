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

wh_is_intercept = find( ~any(diff(X)) ); % which column is the intercept?

switch meth
    
    case 'which'
        X = wh_is_intercept;
        
    case 'remove'
        
        if ~isempty( wh_is_intercept ), X(:, wh_is_intercept) = []; end
        
    case 'add'
        
        X(:, end+1) = 1;
        
    case 'end'
        
        X = intercept(X, 'remove');
        X = intercept(X, 'add');
        
    otherwise
        
        error('Unknown method');
        
end

end
        
