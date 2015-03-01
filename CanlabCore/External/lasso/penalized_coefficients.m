function res = penalized_coefficients(lasso_object, indices, alignment)
% -------------------------------------------------------------------------
% function res = penalized_coefficients(penalized_object, indices, alignment)
% -------------------------------------------------------------------------
% PURPOSE: This function returns the LASSO coefficients for which the 
%          chosen index equals those listed in indices.
%          The alignment index can be chosen to be:
%          lambdas, normalized L1 norms and L1 norms;
% -------------------------------------------------------------------------
% INPUTS:
% lasso_object: a structure as the one returned by the lasso routine:
% indices:      index of the points at which interpolation is wanted
% alignment:    one string of ('normalized_l1', 'l1'm 'lambdas')
%               indicating which index to use for alignment
% -------------------------------------------------------------------------
% OUTPUTS:
% res:          a structure containing  the coefficients and the 
%               regularization parameter at the interpolated points
%               intercept:
%               nbeta:
%               beta:
%               lambda:
% -------------------------------------------------------------------------
% Author: Guilherme V. Rocha
%         Department of Statistics
%         University of California, Berkeley
%         gvrocha@stat.berkeley.edu, gvrocha@gmail.com
% 2006/09
% -------------------------------------------------------------------------
% See also: LASSO, BLASSOL2

% 0. Checking input parameters:
%==========================================================================
if nargin < 3
  alignment = 'normalized_penalty';
end;

% 1. Interpolates according to the chosen L1 norm (normalized x non normalized)
%==========================================================================
switch(lower(alignment))
  case{'normalized_penalty'}
    sizes = lasso_object.npenalty;
  case{'penalty'}
    sizes = lasso_object.penalty;
  case{'lambda'}
    sizes = lasso_object.lambda;
  otherwise
    error('Unrecognized alignment string');
end;

% Takes care of the case when there are entries with repeated L1s:
[unique_sizes, I, J]   = unique(sizes);

% If one of the points to be calculated is beyond the end of the path, use
% the end of path coeffcients for it:
indices       = min(indices, max(sizes));
indices       = reshape(indices, prod(size(indices)), 1);
res.intercept = interp1(unique_sizes, lasso_object.intercept(I), indices);
res.beta      = interp1(unique_sizes, lasso_object.beta(I,:), indices);
res.nbeta     = interp1(unique_sizes, lasso_object.nbeta(I,:), indices);
if(~isempty(lasso_object.lambda))
  res.lambda    = interp1(unique_sizes, lasso_object.lambda(I,:), indices);
end;
res.indices   = indices';
