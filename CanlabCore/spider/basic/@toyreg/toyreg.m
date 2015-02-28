function a = toyreg(hyper) 

%====================================================================
% TOY data for Regression problems
%==================================================================== 
% A=toyreg(H) returns a toy regression object initialized with hyperparameters H. 
%
% This object can generate toy data for regression with
% N input dimensions and O output dimensions. 
% If "nonlinearity" is specified it uses a map object to calculate the
% outputs.
% 
% Hyperparameters, and their defaults
%  l=100          -- examples
%  noiselevel     -- 0
%  n=1            -- input dimension
%  o=1            -- output dimension (possibly overriden by outpmap)
%  nonlinear=0    -- if 1 uses "outpmap" to calculate Y=map(X)
%  outpmap=map;
% 
%  seed=-1       -- random seed used to generate it, if -1 do not
%                   set seed
% 
%
% Methods:
%  generate,train,test

  a.l=100;
  a.noiselevel = 0;
  a.input_noiselevel = 0;
  a.n=1;
  a.o=[];
  a.nonlinear=0;
  a.W=[];
  a.outpmap=map;

  a.seed=-1;
  
  
  p=algorithm('toyreg');
  a= class(a,'toyreg',p);

  if nargin==1
    eval_hyper;
  end  
  
  a.W=randn(a.o,a.n);

