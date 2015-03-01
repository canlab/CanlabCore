clear all;
options.trace = 1;

n_vars = 10;
n_obs  = 100;
beta   = 10.^(-0.3*[1:n_vars])';

rho = 0.85;

correls = rho.^abs([-n_vars:n_vars]);

for i = 1:n_vars
  S(i,:) = correls(n_vars+1-(i-1):2*n_vars-(i-1));
end;
R = chol(S);

x = randn(n_obs, n_vars)*R;
y = x*beta + randn(n_obs, 1);

ttt = lasso(y, x, options);