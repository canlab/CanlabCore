function a = gridsel(alg,hyper)   

%==============================================================      
% GRIDSEL model selection object via grid search   
%==============================================================    
% A=GRIDSEL(A,H) returns a gridsel object initialized with  
% algorithms A and hyperparameters H.   
%  
% Finds the best of the algorithm from the set A, and trains  
% and stores that model.  
%    
% Hyperparameters.  
%  
%  a.child=A                -- methods to evaluate  
%  a.loss='class_loss'      -- loss measure to use  
%  a.score=cv('folds=5')    -- method of evaluating algorithms  
%    
% Model  
%  a.scores=[]      % score of all of methods tried  
%  a.best_score=[]  % score of best methods tried   
%  a.best_index=[]  % index of best methods tried   
%  a.best=[]        % learnt model of best method tried   
%   
% Example:
% % train 3 svms with C=1,2,3 and validate with 3 fold cross validation
% [r,a]=train(gridsel(param(svm,'C',[1,2,3]),{'score=cv;score.folds=3'}),gen(toy)) ;
% 
% Methods:  
%  train, test   
%=============================================================
% Reference : 
% Author    : 
% Link      : 
%=============================================================
   
  
  % hypers   
  a.child=alg;     % original algorithm to start with
  a.score=cv; a.score.folds=5;   
  a.loss='class_loss';  
    
  % model   
  a.all_algs=[];   % all trained models
  a.scores=[];     % score of all of methods tried  
  a.best_score=[]; % score of best method tried   
  a.best_index=[]; % index of best method tried   
  a.best=[];       % learnt model of best method tried (without cv, etc.) 
  
    
  p=algorithm('gridsel');  
  a= class(a,'gridsel',p);  
   
  if nargin==2,  
    eval_hyper;  
  end;  
  
  
