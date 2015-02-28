function [y_pred, y_residuals] = lasso_predict(lasso_obj, x, indices, options)
% -------------------------------------------------------------------------
% function y_pred = lasso_predict(a_lasso_object, x, indices, options)
% -------------------------------------------------------------------------
% PURPOSE: This function returns the estimates for the model contained in
%          blasso_object at the points listed in the matrix x
% -------------------------------------------------------------------------
% INPUTS:
% x:      is a n x p vector where n is the number of points to be predicted
%         and p is the dimension of the regressors
% a_lasso_object: a lasso object as the one returned by the blassol2 or 
%                 lasso functions
%                 it must contain an intercepts and a betas field
% indices: positions along regularization path where predictions will be computed.
%          If left blank, the breakpoints in lasso_obj will be used.
% options: a structure containing:
%          alignment: one of 'penalty', 'normalized_penalty', 'lambda' (default is 'normalized_penalty')
% -------------------------------------------------------------------------
% OUTPUTS:
% y_pred:  the predicted values for the dependent variable along the 
%          regularization path
% -------------------------------------------------------------------------
% Author: Guilherme V. Rocha
%         Department of Statistics
%         University of California, Berkeley
%         gvrocha@stat.berkeley.edu, gvrocha@gmail.com
% 2006/09
% -------------------------------------------------------------------------
% See also: BLASSOL2
% -------------------------------------------------------------------------

if nargin <4;
  options = [];
end;
if nargin < 1;
  error('No parameters defined')
end;
compute_residuals = 0;
if((nargin <2)|(isempty(x)))
  x = repmat(lasso_obj.xmean, lasso_obj.sample_size, 1) + ...
      repmat(lasso_obj.xscale, lasso_obj.sample_size, 1).*lasso_obj.nx;
  y = lasso_obj.ymean + lasso_obj.yscale*lasso_obj.ny;
  compute_residuals = 1;
end;

if(~isfield(options, 'alignment'))
  options.alignment = 'normalized_penalty';
end;

if(~isempty(indices))
  res = lasso_coefficients(lasso_obj, indices, options.alignment); 
else
  res.intercept = lasso_obj.intercept;
  res.beta      = lasso_obj.beta;
end;
n_path = size(res.beta, 1);
n_obs  = size(x, 1);
y_pred = repmat(res.intercept', n_obs, 1) + x*res.beta';
if(compute_residuals)
  y_residuals = y_pred-y;
else
  y_residuals = [];
end;
