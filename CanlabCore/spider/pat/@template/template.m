
function a = template(hyper) 
   
% TEMPLATE template object 
% AUTHOR:  <your name>, <your email>
%  
% A=TEMPLATE(H) returns a template object initialized with hyperparameters H. 
%
%  The template object is a simple linear support vector machine, with
%   soft margin hyperparameter C. Please use this template to
%   implement your own algorithms.
%
% Hyperparameters, and their defaults
%  C=Inf                -- the soft margin C parameter
% 
% Model
%  alpha                -- the weights
%  b0                   -- the threshold
%  Xsv                  -- the Support Vectors
%
% Methods:
%  train, test 
% 
% Example: 
%   [r,a]=train(template,gen(toy))
% ========================================================================
% Refernce  : 
% Author    : 
% Link      : 
% ========================================================================


  a.C=Inf;
  
  % model 
  a.alpha=[];
  a.b0=0;
  a.Xsv=[];
  
  p=algorithm('template');
  a= class(a,'template',p);
 
  if nargin==1,
    eval_hyper;
  end;
