
function a = svr(hyper) 

%=========================================================================   
% SVR svr object
%=========================================================================   
% A=SVR(H) returns a svr object initialized with hyperparameters H. 
%
% 
%  Hyperparameters, and their defaults
%
%   C=Inf;               -- the soft margin C parameter
%   optimizer='andre';  -- other choices={quadprog,andre,libsvm,svmtorch,sparse}
%   alpha_cutoff=-1;     -- keep alphas with abs(a_i)>alpha_cutoff
%   nu = 0;              -- nu parameter of a nu svr (different 
%                           from zero implies the nu-SVR is used
%                           otherwise, the epsilon-SVR is used)
%   child=kernel;        -- the kernel is stored as a member called "child"
%   epsilon=0.1;         -- the value of epsilon in the epsilon 
%                           insensitive loss function during learning
%   use_signed_output=0; -- set to 1 implies that the svr is used
%                           in classification (+1/-1 outputs), set to 0 
%                           implies that the svr is used in regression
%  Model
%
%   alpha                -- the weights
%   b0                   -- the threshold
%   Xsv                  -- the Support Vectors
% 
%  Methods:
%
%   train, test, get_w 
%=========================================================================
% Reference : A tutorial on Support Vector Regression
% Author    : Alex J. Smola , Bernhard Schölkopf
% Link      : http://citeseer.ist.psu.edu/smola98tutorial.html
%=========================================================================

  % hyperparams 
  a.C=Inf;
  
  % NO RIDGE IS USED IN SVR
%  a.ridge=1e-13; 
%  a.balanced_ridge=0;
  a.child=kernel;
  a.optimizer='libsvm';
  a.epsilon = 0.1;
  a.nu = 0;
  a.alpha_cutoff=-1;
 
  % model 
  a.alpha=[];
  a.b0=0;
  a.Xsv=[];
  
  

  
  p=algorithm('svr');
  a= class(a,'svr',p);
  a.algorithm.use_signed_output=0;
 
  if nargin==1,
    eval_hyper;
  end;
  
  
