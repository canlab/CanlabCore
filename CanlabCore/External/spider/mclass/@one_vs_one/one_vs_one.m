function a = one_vs_one(c,hyper) 

%=========================================================================
% ONE_VS_ONE one_vs_one object
%========================================================================= 
% A=ONE_VS_ONE(C,H) returns an one_vs_one object which trains several 
% algorithm C on pairwise problems of class i against class j, 
% and is initialized with hyperparameters H. The classifiers are combined
% by outputting the class with the most votes.
%
% Model
%  child=svm            -- classifier to use for each sub-problem
%
% Methods:
%  train, test, get_w  
%=========================================================================
% Reference : Multi-class Support Vector Machines 
% Author    : Jason Weston , C. Watkins
% Link      : http://citeseer.ist.psu.edu/8884.html
%=========================================================================  
    
  % model
  a.nrofclasses=[];
  if nargin==0
    a.child={svm};  
  else
    a.child=c;   %% algorithms to use  
    if ~isa(c,'cell') 
        a.child={c}; 
    end; 
  end
  
  %% set to use unsigned output if possible 
  for i=1:length(a.child)
    a.child{i}.algorithm.use_signed_output=0;
  end
    
  p=algorithm('one-vs-one');
  a= class(a,'one_vs_one',p);
  %hyperparams are supplied
  if nargin==2
    eval_hyper;    
  end
  
  
