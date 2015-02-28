function a = one_vs_rest(c,hyper) 

%=========================================================================    
% ONE_VS_REST one_vs_rest object
%========================================================================= 
% A=ONE_VS_REST(C,H) returns an one_vs_rest object which trains algorithm 
% C on sub problems of type class i versus all other classes, 
% and is initialized with hyperparameters H. The classifiers are combined
% by outputting the class with the largest positive output (this assumes
% that classifiers output real values indicating confidence rather than
% just +1,-1)
%
% Model
%  child=svm            -- classifier to use for each sub-problem
%
% Methods:
%  train, test, get_w  
% Example:
% 
%c1=[-1,1];c2=[1,1];c3=[0,-1];
% X1= randn(50,2)+repmat(c1,50,1);
% X2= randn(50,2)+repmat(c2,50,1);
% X3= randn(50,2)+repmat(c3,50,1);
% % note the class label format!
% Y1= [ones(50,1),-ones(50,1),-ones(50,1)];
% Y2= [-ones(50,1),ones(50,1),-ones(50,1)];
% Y3= [-ones(50,1),-ones(50,1),ones(50,1)];
% 
% d=data([X1;X2;X3],[Y1;Y2;Y3]);
% 
% [r,a]=train(one_vs_rest(svm(kernel('rbf',2))),d)
% % Test class centers
% dtest=data([c1;c2;c3]);
% rtest=test(a,dtest)
%=========================================================================
% Reference : Multi-class Support Vector Machines 
% Author    : Jason Weston , C. Watkins
% Link      : http://citeseer.ist.psu.edu/8884.html
%=========================================================================
    
  % model
  if nargin==0
    a.child={svm};  
  else
    a.child=c;   %% algorithms to use  
    if ~isa(c,'cell') a.child={c}; end; 
  end
  
  %% set to use unsigned output if possible 
  for i=1:length(a.child)
    a.child{i}.algorithm.use_signed_output=0;
  end
    
  p=algorithm('one-vs-rest');
  a= class(a,'one_vs_rest',p);

  if nargin==2
    eval_hyper;    %evaluate hyperparams
  end
  
  
