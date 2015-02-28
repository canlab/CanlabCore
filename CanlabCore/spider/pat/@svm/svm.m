function a = svm(hyper) 
%=============================================================================
% SVM Support Vector Machine object             
%=============================================================================  
% a=svm(hyperParam) 
%
% Generates a svm object with given hyperparameters.
%
%
%   Hyperparameters (with defaults)
%   child=kernel         -- the kernel is stored as a member called "child"
%   C=Inf                -- the soft margin C parameter
%   ridge=1e-13          -- a ridge on the kernel
%   balanced_ridge=0     -- for unbalanced data
%   nu = 0               -- Schoelkopf's nu svm parameter
%   optimizer='default'  -- other choices={andre,quadprog,svmlight,
%                                          libsvm,svmtorch(linux only)}
%                           For "libsvm" optimizer you can specify the used cache size
%                           by the global variable "libsvm_cachesize".
%   alpha_cutoff=-1;     -- keep alphas with abs(a_i)>alpha_cutoff
%                           default keeps all alphas, another
%                           reasonable choice is e.g alpha_cutoff=1e-5 to remove
%                           zero alphas (i.e non-SVs) to speed up computations.
% 
%   Model
%    alpha               -- the weights
%    b0                  -- the threshold
%    Xsv                 -- the Support Vectors
%
% Methods:
%  train, test, get_w 
%
% Example:
%
%  d=gen(spiral({'m=200','n=2','noise=0.35'}));
%  [r,a]=train(cv(svm({kernel('rbf',1),'optimizer="andre"'})),d)
%  plot(a{1})
%
%=============================================================================
% Reference : A Tutorial on Support Vector Machines for Pattern Recognition  
% Author    : Christopher J. C. Burges
% Link      : http://citeseer.ist.psu.edu/burges98tutorial.html
%=============================================================================

  %<<------hyperparam initialisation------------->> 
  a.child=kernel;
  a.C=Inf;
  a.ridge=1e-13;  
  a.balanced_ridge=0;
  a.nu = 0;
  a.optimizer='default';
  a.alpha_cutoff=-1;
  
  
  % <<-------------model----------------->> 
  a.alpha=[];
  a.b0=0;
  a.Xsv=[];
  a.nob=0;
  
  algoType=algorithm('svm');
  a= class(a,'svm',algoType);

  a.algorithm.alias={'kern','child'}; % kernel aliases
  
 if nargin==1,
    eval_hyper;
 end;





