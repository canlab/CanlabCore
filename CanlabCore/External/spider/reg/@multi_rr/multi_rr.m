  
function a = multi_rr(hyper)   
 
 %========================================================================     
% MULTI_RR (possibly multiple output) ridge regression object  
%=========================================================================    
% A=MULTI_RR(H) returns an object initialized with hyperparameters H.   
%  
% Performs ridge regression independently for each column we wish  
% to learn, which only means inverting a single (kernel) matrix  
%  
% Hyperparameters, and their defaults  
%  ridge=1e-13;   -- a ridge on the kernel  
%  child=kernel;  -- the kernel is stored as a member called "child"  
%  use_kernels=1; -- if this is set to 0 a linear ridge regression or 
%                    an empirical kernel map (useful for reduced sets
%                    of centers)
%  indices = []   -- indices of a reduced set of centers to be used  
%                       for learning. ([] means use all training set)  
%  use_b=1        -- find a threshold, otherwise fix to 0   
%  
% Model  
%  alpha          -- the weights  
%  Xsv            -- centers 
%  
% Methods:  
%  train, test, get_w   
%  loo : calculate leave one out predictions (with empirical kernel, or linear)
%=========================================================================    
% Reference : Ridge Regression Learning Algorithm in Dual Variables
% Author    : C. Saunders , A. Gammerman and V. Vovk
% Link      : http://citeseer.ist.psu.edu/saunders98ridge.html
%=========================================================================  
  
  %hyperparams   
  a.ridge=10^(-13);   
  a.child = kernel;  
  a.indices = []; 
  a.use_kernels=1;
  a.use_b=1;  
     
  % model   
  a.alpha = [];  
  a.b0=0;  
  a.Xsv=[];  
    
  p=algorithm('multi_rr');  
  p.use_signed_output=0;  
  a= class(a,'multi_rr',p);  
   
  if nargin==1,  
    eval_hyper;  
  end;  
  
    
  
