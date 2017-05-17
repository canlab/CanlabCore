  
function a = kpca(hyper)   

%================================================================================   
% KPCA kpca object - Kernel Principal Components Analysis  
%================================================================================    
% A=KPCA(H) returns a kpca object initialized with hyperparameters H.   
%  
% Hyperparameters, and their defaults  
%  feat=0;              -- number of features, default 0 means all via rank(K) 
%  center_data=1;       -- if data is to be centered in feature space
%  child=linear         -- child stores the kernel. Default is the linear
%                          kernel and therefore normal pca. 
%                          NOTE: This has changed. The old version was
%                          assuming a kernel matrix as data. In order to
%                          simulate the old behaviour use custom kernel.
% Model  
%  e_val                -- the eigenvectors  
%  e_vec                -- the eigenvalues  
%  dat                  -- training data (that we extracted from)  
%  
% Methods:  
%  train, test  
%================================================================================
% Reference : Nonlinear component analysis as a kernel eigenvalue problem
% Author    : B. Schölkopf, A. Smola, and K.-R. Müller
% Link      : http://www.kernel-machines.org/papers/nlpca.ps.gz
%================================================================================
    
  %hyperparams   
  a.feat=0;  
  a.center_data = 1;
  a.child=kernel('linear');  
    
  % model   
  a.e_vec=[]; % eigenvectors  
  a.e_val=0;  % eigenvalues  
  a.dat=[];  
  a.Kt=[];
  p=algorithm('kpca');  
  a= class(a,'kpca',p);  
   
  if nargin==1,  
    eval_hyper;  
  end;  
 
