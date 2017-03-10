function a = gproc(hyper) 
   
%=============================================================================
% GPROC Gaussian-Process object
%=============================================================================  
%   
% A=GPROC(hyper) returns a GPROC object initialized with hyperparameters hyper. 
% minimize log-likelihood for length runs
%
% 
%  Hyperparameters, and their defaults
%
% 
%   child=kernel;        -- the kernel is stored as a member called "child"
%   length=50            -- default maximum iteration of line searches.
%  
%  Model
%  H			-- covariance hyperparameters
% 
%  Methods:
%  d=gen(toyreg({'o=1','n=2','l=200'}))
%  
%  [r,a]=train(gproc,get(d,1:100))
%  r=test(a,get(d,101:200))
%  % Note that r.X has 2 an extra output columns! 
%
%   train, test
%=============================================================================
% Reference : Introduction to gaussian processes
% Author    : David J.C. MacKay
% Link      : http://www.inference.phy.cam.ac.uk/mackay/gpB.ps.gz
%=============================================================================

  % hyperparams 
  a.input= [];
  a.target = [];
  a.length=50;
  a.H=[];
  a.pred_var=0;
  % model     
    
  p=algorithm('gproc');
  a= class(a,'gproc',p);
  a.algorithm.use_signed_output=0;
 
  if nargin==1
    eval_hyper;
  end;
  
  
