
function a = gauss(i1,i2,i3) 

%====================================================================   
% GAUSS gauss/gaussian density estimation object
%==================================================================== 
% A=GAUSS([M],[S],H) returns a gauss object initialized with mean M,
% std S and hyperparameters H. Can also be called with GAUSS(H)
%
% Training will try to fit a single Gaussian to the data.
% Testing will return the probability estimate, or if passed
% an empty dataset will generate new data according to the 
% density learnt.
%
% Hyperparameters:
%  l=50       -- number of data points to generate if asked to generate   
%  assume=[]  -- can make assumptions: {'diag_cov','equal_cov'}
%                where the matrix is diagonal ('diag_cov'), 
%                or diagonal with all elements equal ('equal_cov')
%  
% Model:
%  mean=0  -- centre of gauss distbn
%  cov=1   -- cov matrix, if single value assumes diagonal
%             of matrix all has same values, if vector assumes it
%             is the diagonal of the matrix, with 0s everywhere else
%
% Methods:
%  train, test, generate
%
% Examples:
%  gen(gauss)                                   % generate data
%  gen(gauss('l=5;mean=[1 1];cov=0.01;'))       % generate data
%
%  d=test(gauss([1 5],[1 -0.4 ; -0.4 1],'l=300')); 
%  [d2 a]=train(gauss('l=200'),d);
%  d2=test(a); %% train a gauss on some (gauss) data and then
%              %% try to generate similar data
%  hold off; plot(d.X(:,1),d.X(:,2),'o')
%  hold on;  plot(d2.X(:,1),d2.X(:,2),'rx')
%====================================================================
% Reference : chapter 2 (Richard O. Duda and Peter E. Hart) Bayesian Decision Theory
% Author    : Richard O. Duda , Peter E. Hart
% Link      : http://www.amazon.com/exec/obidos/tg/detail/-/0471056693/002-6279399-2828812?v=glance
%====================================================================  
  
  a.l=50;
  a.mean=0;
  a.cov=1;
  a.assume=[];
  
  p=algorithm('gauss');
  a= class(a,'gauss',p);
 
  
  if nargin==1 
    if ischar(i1)
      hyper=i1; eval_hyper; return;
    else
      a.mean=i1; 
    end
  end;
  
  if nargin==2 
    a.mean=i1; 
    if ischar(i2)
      hyper=i2; eval_hyper; return;
    else
      a.cov=i2;
    end
  end;
  
  if nargin==3 
    a.mean=i1;
    a.cov=i2;
    hyper=i3; eval_hyper; 
  end;
  






