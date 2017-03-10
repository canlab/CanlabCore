
function a = cv(algo,hyper) 

%=====================================================================================
%  Cross validation object      
%=====================================================================================   
%  a=cv(algo,hyperParam) Returns a cv object on algorithm algo using with 
%                        given Hyperparameters. 
% 
%  Possible hyperparameters (with defaults):
%   folds=5               - number of folds
%   repeats=1             - can do n*CV for reduced variance
%   balanced=0            - determines if cv shall be balanced  (same number of
%                            positives in each fold)
%   store_all=1           - determines if models trained in all folds shall be stored   
%   output_train_error=0  - determines wheter to output training error on each fold
%                            (else cv error, i.e  test left out fold as default)
%   train_on_fold=0       - determines to  test on left out fold and train on the
%                            rest (set to true for the opposite).
%   store_trialbytrial=0  - whether to store the field a.trialbytrial, the rows of
%                            which are [i, foldnumber, y_i, f(x_i)] for each index
%                            i of a data point on which the algorithm was tested
%
%  Model:
%   child                 - stored in child algorithm 
%
%  Methods:
%   train                 - selfexplanatory
%   test                  - selfexplanatory
%=====================================================================================
% Reference : 
% Author    : 
% Link      : 
%=====================================================================================
  
  % <<---- Initialisation of Hyperparams ---->>
  a.folds=5;
  a.repeats=1;            
  a.balanced=0; 
  a.store_all=1;    
  a.output_train_error=0; 
  a.train_on_fold=0;      
  a.store_trialbytrial=0;
  a.trialbytrial=[];

  % <<---- Initialisation of model ---->>
  if nargin>0
   a.child{1}=algo;
  else
   a.child{1}=svm; 
  end 
    
  % <<---- convert {} to group({}) ---->>
  if isa(a.child{1},'cell') 
      a.child{1}=group(a.child{1}); 
  end;
  
  
  p=algorithm('cv');
  a= class(a,'cv',p);
  if nargin==2,
    eval_hyper;
  end;
  
  if a.store_all
    for i=2:a.folds*a.repeats
     a.child{i}=a.child{1}; 
    end
  end
   
