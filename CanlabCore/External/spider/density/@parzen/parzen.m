
function a = parzen(hyper) 

%=========================================================================   
% PARZEN parzen's windows density estimation object
%========================================================================= 
% A=PARZEN(H) returns a parzen object initialized with hyperparameters H. 
%
% Kernel density estimation.
%
% Model
%  child=kernel('Gaussian') -- the kernel to use
%
% Methods:
%  train, test 
%
% Sample demo code:
% -----------------
%d=data(randn(30,1));
%dtst=data([-5:0.1:5]');
%figure(1);clf;
%for j=1:9
% subplot(3,3,j);
% plot(d.X,d.X*0,'x'); 
% a=parzen; a.kerparam=1/(2^(j-2)); 
% [r a]=train(a,d); r=test(a,dtst);
% hold on ; [u ind]=sort(dtst.X); 
% plot(dtst.X(ind,:),r.X(ind,:),'-');
% title(['\sigma=1/2^{-' num2str(j) '}'])
%end
%=========================================================================
% Reference : chapter 4 (Richard O. Duda and Peter E. Hart) Parzen Windows
% Author    : Richard O. Duda , Peter E. Hart
% Link      : http://www.amazon.com/exec/obidos/tg/detail/-/0471056693/002-6279399-2828812?v=glance
%=========================================================================

  % model 
  a.child=kernel('Gaussian',10);
  
  p=algorithm('parzen');
  a= class(a,'parzen',p);
 
  if nargin==1,
    eval_hyper;
  end;
