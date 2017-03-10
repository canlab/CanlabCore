
function a = toy(hyper) 

%====================================================================
% TOY toy data generation object
%==================================================================== 
% A=TOY(H) returns a toy object initialized with hyperparameters H. 
%
% This generates toy data from a uniform distribution in (n)
% dimensions, only (relevant_n) of which are used to label the data
% (the last few features n-relevant_n:n are relevant).
% 
% Hyperparameters, and their defaults
%  l=30          -- examples
%  n=30          -- dimensions
%  relevant_n=10 -- number of dimensions relevant to output
%  seed=-1       -- random seed used to generate it, if -1 do not
%                   set seed
% 
% Model
%  w             -- stores label model
%
% Methods:
%  generate,train,test

  a.l=50;
  a.n=30;
  a.relevant_n=10;
  a.seed=-1;
  
%  a.uniform_x=1;
%  a.random_weights=1; 
%  a.weights=[];
  
  p=algorithm('toy');
  a= class(a,'toy',p);
%  a.weights=ones(1,a.n); 

  if nargin==1
    eval_hyper;
  end  
