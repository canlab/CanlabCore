function a = bayessel(A,hyper)

%==============================================================================================================
% Bayesian model selection for SVM/SVR following Kwok et al. , 
% Proceedings of the European Symposium on Artificial Neural Networks (ESANN), pp.177-182, Bruges, Belgium, April 1999.
% for selecting the parameter C and ranking models
%==============================================================================================================
%
% inputs:
%   A		    - SVM/SVR
%
% hyperparameters and their defaults:
%   use_balanced_C=0  - adapt balanced C for SVM
%   type='L1'	      - use l1-SVM
%
% Methods:
%  train    - optimizes soft margin parameter for given data
%  test     - 
%
% outputs:
%   pbest      - best parameter C
%   posterior  - posterior probability for the model (given the data)
% Example:
% 
%  [r,a]=train(bayessel(svm),gen(toy))
%  [a.pbest,a.posterior]
%==============================================================================================================
% Reference : (The Evidence Framework Applied to SVM / Bayesian SVR)
% Author    :  Kwok et al
% Link      : Proceedings of the European Symposium on Artificial Neural Networks (ESANN), pp.177-182
%==============================================================================================================  
  
  %hyperparams 
  a.child=A;
  a.use_balanced_C=0;  
  a.type='L1';

  % model   
  a.pbest = [];
  a.posterior=-inf;
  
  p=algorithm('bayessel');
  a= class(a,'bayessel',p);
 

  if nargin==2,
    eval_hyper;
  end;
