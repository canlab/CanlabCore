function a = kvq(hyper)   
%================================================================================   
% KVQ kvq object - Kernel Vector Quantization
%================================================================================    
% A=KVQ(H) returns a kvq object initialized with hyperparameters H.   
%  
% The kvq object computes a learning vector quantization with guaranteed distortion 
% bounds. It does so by minimizing the L1 norm (parameter optmethod='l1lvq' / c.f. cited
% reference below) or an approximation of the L0 norm (parameter optmethod='l0arom').
% 
% It also implements two variants of LVQ for labeled data: discriminative and shared
%  LVQ (set parameter mode to 'discriminative' or 'shared' respectively). In 
% discriminative LVQ a vector can only become a codebook vector if it has at least 
% distance delta (parameter distd) to any examples of the other class. In mode shared an 
% example can only become a codebook vector if it has at most distance delta to a example
% of the opposite class. 
%
% Hyperparameters, and their defaults  
%  dist=1;          -- allowed point to point distortion.
%  distd=3			-- discriminative or shared radius
%  child=distance         -- child stores the distance object.
%  a.cutoff = 1./3.       -- cutoff value for importance values
%  a.return_indices = 0   -- return indices of points instead of samples
%  a.tol = 1e-6			  -- tolerance parameter for matlab's LINPROG linear program optimizer
%  a.test_on_trainingset = 1 -- test the algorithm on the training data
%
% Model  
% a.keep                -- kept data points -- store data points kept for model       
% a.alpha               -- importance factors
% 
%
% Methods:  
%  train, test,circplot
%
% Example: 
%  d=gen(toy('n=2'));
%  [r,a]=train(kvq,d)
%  plot(a,d)
%================================================================================
% Reference : A kernel approach vor vector quantization with guaranteed distortion bounds
% Author    : M. Tipping, B. Schölkopf
% Link      : ftp://ftp.research.microsoft.com/users/mtipping/aistats01.ps.gz
%================================================================================
    
  %hyperparams   
  a.dist=1;  % epsilon
  a.distd = 3; % delta
  a.child=distance;  

  % cutoff value for samples.
  a.cutoff=1./3.;
    
  a.mode = 'standard';
  
  % return indices of points instead of samples.
  a.return_indices=0;
  % model   
  a.alpha =[];
  a.keep=[];

  a.optimizer = 'linprog';
  a.optmethod= 'l0arom';
  a.tolfun = 1e-8;
  a.tol = 1e-6;
  a.test_on_trainingset = 1;
  
  p=algorithm('kvq');  
  a= class(a,'kvq',p);  
  
   
  if nargin==1,  
    eval_hyper;  
  end;  
 
