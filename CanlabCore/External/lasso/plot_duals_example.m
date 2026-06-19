clear all;
options.trace = 1;

% n_vars: Number of variables.
% n_obs: Number of observations.
% beta: Coefficient vector for the variables.
% rho: Correlation coefficient for the variables.
n_vars = 10;
n_obs  = 100;
beta   = 10.^(-0.3*[1:n_vars])';
rho = 0.85;

% correls: Creates a vector of correlations.
correls = rho.^abs([-n_vars:n_vars]);

% S: Initializes the correlation matrix.
% The loop fills the correlation matrix S using the correlations vector.
for i = 1:n_vars
  S(i,:) = correls(n_vars+1-(i-1):2*n_vars-(i-1));
end;

% Computes the Cholesky decomposition of the correlation matrix.
% The Cholesky decomposition is a mathematical technique used to decompose 
% a symmetric positive-definite matrix into the product of a lower triangular matrix and its transpose.
% resulting in a lower triangular matrix R.
R = chol(S);

% Multiplying by R transforms these independent random variables into a set 
% of correlated random variables x that have the correlation structure defined by S.
x = randn(n_obs, n_vars) * R: % Generates the correlated predictor variables.
y = x * beta + randn(n_obs, 1): % Generates the response variable. The correction here is using x * beta to ensure y is a vector.

% Performs LASSO regression using the predictor matrix x and response vector y.
% ttt = lasso(y, x, options); This is a bug. X and Y were reversed. Needs
% 'Options', prior to options.
ttt = lasso(x, y, 'Options', options);