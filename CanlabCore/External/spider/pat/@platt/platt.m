function a = platt(c,hyper)
%============================================================================ 
% Probabilistic outputs for SVM
%============================================================================
% A=PLATT(C,H) returns a platt object initialized with hyperparameters H and
%             based on algorithm C. 
% Converts a real valued margin producing pattern recognition
% algorithm into a conditional probability estimator. It thus
% requires a pattern recognition approach that can produce a real
% valued output before thresholding, e.g SVMs. It achieves this
% by finding the coefficients of a sigmoid by using cross validation.
%
% Hyperparameters, and their defaults
%
%  child = svm          -- real valued margin algorithm to be used 
%  folds = 3            -- number of cv folds
%
% Model:
%  A                -- 
%  B                -- Parameters of the sigmoid: 1/(1+exp(A*f(x)+B)
%  
% Methods:
%  training, testing 
% 
% Example:
% [rr,a]=train(platt(svm),gen(toy))
% % rr.X  contains now the probability of being class 2
% % predictions 
% [sign(rr.X-0.5),rr.Y]
%============================================================================
% Reference : Probabilistic Outputs for Support Vector Likelihood Methods Machines and Comparisons to Regularized
% Author    : John Platt
% Link      : http://research.microsoft.com/users/jplatt/SVMprob.ps.gz
%============================================================================
  
  % model 
  
  a.child = svm;
  a.folds = 3;
  a.B=[];
  a.A=[];
  
  if nargin==0
    a.child=svm;  
  else 
    a.child=c; %% algorithm to use  
  end
  
  
  p=algorithm('platt');
  p.use_signed_output=0;
  a= class(a,'platt',p);
 

  if nargin==2,
    eval_hyper;
  end;







