
function a = bayes(i1,i2) 

%===============================================================  
% BAYES simple bayesian classification object
%===============================================================
% A=BAYES(C,H) returns a bayes object initialized with density
% estimator C and hyperparameters H. 
%
% Training will try to fit the estimator C to each class of 
% pattern recognition data. C can be an array of algorithms, one
% for each class.
% Testing will return the class with the highest probability
% estimate, or if passed the empty dataset will generate new 
% class data according to the densities learnt 
%
% Hyperparameters:
%  l=50            -- number of data points to generate if asked to generate   
%  train_priors=1  -- whether to adjust priors or not according to data
%  
% Model:
%  child=gauss     -- array of underlying density estimators
%  prior=1         -- prior probabilities for each class, default equal
%
% Methods: 
%  train, test, generate
%
% Examples:
%   get_mean(train(cv({bayes svm}),toy))
%   get_mean(train(cv({bayes(gauss('assume="diag_cov"')) svm}),toy))
%   gen(bayes({gauss([-1]) gauss([1])}))
%   d=gen(bayes({gauss([-1 3]) gauss([0 4]) gauss([1 2])}))
%===============================================================
% Reference : chapter 2 (Richard O. Duda and Peter E. Hart) Bayesian Decision Theory
% Author    : Richard O. Duda , Peter E. Hart
% Link      : http://www.amazon.com/exec/obidos/tg/detail/-/0471056693/002-6279399-2828812?v=glance
%===============================================================

  a.l=50;
  a.child=gauss;
  a.prior=1;
  a.train_priors=1; 
  
  p=algorithm('bayes');
  a= class(a,'bayes',p);
 
  if nargin==1 
    if ischar(i1)
      hyper=i1; eval_hyper; 
    else
      a.child=i1; 
    end
  end;
  
  if nargin==2 
    a.child=i1; 
    hyper=i2; eval_hyper; 
  end;
  
  a.prior=ones(1,length(a.child))/length(a.child);





