function a = multi_reg(c,hyper) 

%=========================================================================   
% MULTI_REG object - for regression across multiple outputs 
%=========================================================================  
% A=MULTI_REG(C,H) returns an multi_reg object which uses C as a
% base classifier for each output, and is initialized with hyperparameters H. 
%
% Hyperparameters, and their defaults
%  child=C              -- the base classifier for each output
% 
% Model
%  [stored in children]
%
% Methods:
%  train, test
%
% Example:
%   a=multi_reg(svr('C=10'));
%   d=toyreg('n=10;o=5;');
%   [r a]=train(a,d) 
%
%=========================================================================
% Reference : 
% Author    : 
% Link      : 
%=========================================================================

  a.Q=0; % number of outputs
  % model
  if nargin==0
    a.child={svr};  
  else
    a.child=c;   %% algorithms to use  
    if ~isa(c,'cell') a.child={c}; end; 
  end
  
  p=algorithm('multi_reg');
  a= class(a,'multi_reg',p);

  %hyperparams
  if nargin==2
    eval_hyper;    
  end


