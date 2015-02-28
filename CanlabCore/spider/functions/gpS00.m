function [out1, out2] = gpS00(X, P);
%%function [out1, out2] = gpS00(X, {input, target, test});

% gpS00: Gaussian process regression with "squared negative exponential"
% covariance function and independent Gaussian noise model. Two modes are
% possible: training and prediction: if no test data are given, the function
% returns minus the log likelihood and its partial derivatives with respect to
% the hyperparameters; this mode is used to fit the hyperparameters. If test
% data are given, then (marginal) Gaussian predictions are computed, whose mean
% and (noise free) variance are returned.
%
% usage: [fX dfX] = gpS00(X, input, target)
%    or: [mu S2]  = gpS00(X, input, target, test)
%
% where:
%
%   X      is a (column) vector (of size D+2) of hyperparameters
%   input  is a n by D matrix of training inputs
%   target is a (column) vector (of size n) of targets
%   test   is a nn by D matrix of test inputs
%   fX     is the returned value of minus log likelihood
%   dfX    is a (column) vector (of size D+2) of partial derivatives
%            of minus the log likelihood wrt each of the hyperparameters
%   mu     is a (column) vector (of size nn) of prediced means
%   S2     is a (column) vector (of size nn) of predicted variances
%
% where D is the dimension of the input. The form of the covariance function is
%
% C(x^p,x^q) = v^2 * exp[-(x^p - x^q)'*inv(P)*(x^p - x^q)/2]
%            + u^2 * delta_{p,q}
%
% where the first term is the squared negative exponential and the second term
% with the kronecker delta is the noise contribution. The P matrix is diagonal
% with "Automatic Relevance Determination" (ARD) or "input length scale"
% parameters w_1^2,...,w_D^2; The hyperparameter v is the "signal std dev" and
% u is the "noise std dev". All hyperparameters are collected in the vector X
% as follows:
%
% X = [ log(w_1)
%       log(w_2) 
%        .
%       log(w_D)
%       log(v)
%       log(u) ]
%
% Note: the reason why the log of the parameters are used in X is that this
% often leads to a better conditioned (and unconstrained) optimization problem
% than using the raw hyperparameters themselves.
%
% This function can conveniently be used with the "minimize" function to train
% a Gaussian process:
%
% [X, fX, i] = minimize(X, 'gpS00', length, input, target)
%
% See also: minimize, hybrid
%      
% (C) Copyright 1999 - 2003, Carl Edward Rasmussen (2003-08-08).

input=P{1};
target=P{2};
if length(P) == 3
	test=P{3};
	testmode=1;
else
	testmode=0;
end

[n, D] = size(input);         % number of examples and dimension of input space
input = input ./ repmat(exp(X(1:D))',n,1);

% first, we write out the covariance matrix Q

warning off;

Q = zeros(n,n);
for d = 1:D
  Q = Q + (repmat(input(:,d),1,n)-repmat(input(:,d)',n,1)).^2;
end
Q = exp(2*X(D+1))*exp(-0.5*Q);

if ~testmode   % if no test cases, we compute the negative log likelihood ...

  W = (Q+exp(2*X(D+2))*eye(n))\eye(n);               % W is inv (Q plus noise term)
  invQt = W*target;                               % don't compute determinant..
  [C,flag]=chol(Q+exp(2*X(D+2))*eye(n));
  logdetQ = 2*sum(log(diag(C)));        % ..directly
  out1 = 0.5*logdetQ + 0.5*target'*invQt + 0.5*n*log(2*pi);

  % ... and its partial derivatives

  out2 = zeros(D+2,1);                  % set the size of the derivative vector
  W = W-invQt*invQt';
  Q = W.*Q;
  for d = 1:D
    out2(d) = ...
            sum(sum(Q.*(repmat(input(:,d),1,n)-repmat(input(:,d)',n,1)).^2))/2;
  end 
  out2(D+1) = sum(sum(Q));
  out2(D+2) = trace(W)*exp(2*X(D+2));

else                    % ... otherwise compute (marginal) test predictions ...

  [nn, D] = size(test);     % number of test cases and dimension of input space
  test = test ./ repmat(exp(X(1:D))',nn,1);

  a = zeros(n, nn);    % compute the covariance between training and test cases
  for d = 1:D
    a = a + (repmat(input(:,d),1,nn)-repmat(test(:,d)',n,1)).^2;
  end
  a = exp(2*X(D+1))*exp(-0.5*a);

  % ... write out the desired terms

  if nargout == 1
    out1 = a'*((Q+exp(2*X(D+2))*eye(n))\target);              % predicted means
  else
    invQ = inv(Q+exp(2*X(D+2))*eye(n));
    out1 = a'*(invQ*target);                                  % predicted means
    out2 = exp(2*X(D+1)) - sum(a.*(invQ*a),1)'; % predicted noise-free variance
  end

end
