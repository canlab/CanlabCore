function result = lasso_cv(y, x, k, options)
% -------------------------------------------------------------------------
% function result = lasso_cv(y, x, k, options)
% -------------------------------------------------------------------------
% PURPOSE: fits the parameters of a linear model by using the lasso and 
%          k-folds cross validation
% -------------------------------------------------------------------------
% INPUTS:
% y:         the values of the dependent variable
% x:         the values of the explanatory variables - constants will be 
%            ignored
% k:         determine number of cross validation folds
% n_div: 
% alignment: 
% cv_assignment: allows that the division within groups to be defined by
%                the user. cv_assignment is a nx1 vector where n is the
%                number of available observations. The second columns
%                defines which observations belongs to each group. 
%                If this
%                parameter is used, k is ignored and set to be equal to the
%                number of unique elements of cv_assignment;
% -------------------------------------------------------------------------
% OUTPUTS:
% a structure containing:
% all the fields contained in the strucutre returned by the LASSO function
% (see the lasso function help)  plus the following fields:
% result.CV.intercept    : intercept for the model chosen by CV
% result.CV.betas        : coefficients for the model selected by CV
% result.CV.fit          : Fitted values for the cross-validated s
% result.CV.residuals    : Residuals for the cross-validated s
% result.CV.SSR          : Sum of squared residuals for the cross-validated
%                          estimate
% result.CV.MSE          : minimum cross-validated MSE
% result.CV.cv_assignment: a matrix determining which observations where 
%                          assigned to which CV groups
% result.CV.k            : number of CV folds
% -------------------------------------------------------------------------
% WARNING:
% unless cv_assignment is provided, this function invokes ASSIGN_CV.
% ASSIGN_CV calls randperm affecting the state of the random number gen.
% -------------------------------------------------------------------------
% Author: Guilherme V. Rocha
%         Department of Statistics
%         University of California, Berkeley
%         gvrocha@stat.berkeley.edu, gvrocha@gmail.com
% 2006/09
% -------------------------------------------------------------------------
% See also: LASSO, ASSIGN_CV, RANDPERM

if(nargin < 4)
  options = [];
end;

if(~isfield(options, 'alignment'))
  options.alignment = 'normalized_penalty';
end;
if(~isfield(options, 'cv_assignment'))
  assign_cv_flag = 1;
else
  assign_cv_flag = 0;
  cv_assignment  = options.cv_assignment;
  k              = length(unique(cv_assignment));
end;
if(~isfield(options, 'num_div'))
  options.num_div =100;
end;  
n_div = options.num_div;
alignment = options.alignment;

n_obs           = size(y, 1);
max_index       = -Inf;
if(assign_cv_flag)
  cv_assignment = assign_cv(n_obs, k);
end;

for i = 1:k
  fit_indexes   = find(cv_assignment~=i);
  y_fit         = y(fit_indexes);
  x_fit         = x(fit_indexes,:);
  res(i)        = lasso_rocha(y_fit, x_fit, options);
  switch lower(alignment)
    case 'normalized_penalty'
      max_index = max(max_index, max(res(i).npenalty));
    case 'penalty'
      max_index = max(max_index, max(res(i).penalty));
    case 'lambda'
      max_index = max(max_index, max(res(i).lambda));
  end;
end;  

index_slices = 0:max_index/n_div:max_index;

for i = 1:k
  interpol    = lasso_coefficients(res(i), index_slices, alignment);
  cv_indexes  = find(cv_assignment==i);
  y_cv        = y(cv_indexes);
  x_cv        = x(cv_indexes,:);

  y_pred      = lasso_predict(res(i), x_cv, index_slices, options);
  residuals   = repmat(y_cv, 1, size(y_pred, 2)) - y_pred;
  if(size(residuals, 1)>1)
    MSE(i,:)    = mean(residuals.^2);
  else
    MSE(i,:)    = residuals.^2;
  end;  
end;

[min_MSE, best_cv_MSE] = min(mean(MSE));
best_index             = index_slices(best_cv_MSE);
result                 = lasso_rocha(y, x, options);

result.CV                = lasso_coefficients(result, best_index, alignment);
result.CV.cv_index       = best_cv_MSE;
result.CV.cv_assignment  = cv_assignment;
result.CV.k              = k;
result.CV.min_MSE        = min_MSE;
result.CV.MSEs           = MSE;
result.CV.path           = lasso_coefficients(result, index_slices, alignment);
result.options           = options;
rm_fields                = {'ny', 'nx', 'xtx', 'xty', 'beta', 'intercept', 'penalty', 'xtrs'};
result.CV.fold_results   = rmfield(res, intersect(fieldnames(res), rm_fields));
