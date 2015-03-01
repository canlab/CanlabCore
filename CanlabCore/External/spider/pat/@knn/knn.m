function a = knn(hyper) 
%========================================================================
% KNN k-nearest neighbours object
%========================================================================  
% a = knn(hyperParam)
% 
% Generates a knn object with given hyperparamters.
%
% Hyperparameters (with defaults)
%  k=1                  -- number of neighbours
%  batch=1              -- true if to be computed in batch 
%                           (requires more memory)
%  child=kernel;        -- kernel function (distances computation/ 
%                           default is linear kernel)
%  output_preimage=0    -- whether output index from training sample 
%                           of preimage 
%                          instead of actual label
%
% Model
%  dat                  -- data used for neighbours computation
%
% Methods:
%  train, test
% c1=[2,0];
% c2=[-2,0];
% X1=randn(50,2)+repmat(c1,50,1);
% X2=randn(50,2)+repmat(c2,50,1);
% d=data([X1;X2],[ones(50,1);-ones(50,1)]);
% [r,a]=train(cv(knn({kernel('poly',3),'k=3'})),d); % use knn with 3 number of neighbours and polynomial kernel order 3
% loss(r)
%========================================================================
% Reference : chapter 4 (Richard O. Duda and Peter E. Hart) k-nearest neighbor
% Author    : Richard O. Duda , Peter E. Hart
% Link      : http://www.amazon.com/exec/obidos/tg/detail/-/0471056693/002-6279399-2828812?v=glance
%========================================================================
  
  % <<------------hyperparams initialisation -------->>
  a.k=1;
  a.child=kernel;  
  a.no_train=0;
  a.batch=1;
  a.output_preimage=0;
  
  % <<-------------model-------------------------->> 
  a.dat=[];
    
  algoType=algorithm('knn');
  a= class(a,'knn',algoType);
  if nargin==1,
    eval_hyper;
  end;  
  
  a.algorithm.use_signed_output=0;
   


