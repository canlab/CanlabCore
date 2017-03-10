
function a = one_class_svm(hyper) 

%=============================================================================
% ONE_CLASS_SVM svm object
%============================================================================= 
% A=ONE_CLASS_SVM(H) returns an svm object initialized with hyperparameters H. 
%
%  Learns one class problems, e.g for novelty detection by learnign
%  the support of a density.
%
% Hyperparameters, and their defaults
%  C=Inf                -- the soft margin C parameter
%  optimizer='default'  -- other choices={andre,quadprog,svmlight,libsvm}
%  nu = 0               -- bernhard's nu svm parameter
%  child=kernel         -- the kernel is stored as a member called "child"
% 
% Model
%  alpha                -- the weights
%  b0                   -- the threshold
%  Xsv                  -- the Support Vectors
%
% Methods:
%  train, test, get_w 
%  Example: 
%    % Assume only positive examples are available
%   d=gen(toy('l=100'));
%   d.X=d.X(d.Y==1,:); d.Y=d.Y(d.Y==1,:); % only have one class available
%   [r a]=train(param(one_class_svm('optimizer="quadprog"'),'nu', 2.^[-5:0.5:-1]),d); 
%   loss(test(a,gen(toy)))
%   pause
%   d0=gen(toy2d('2circles','l=200'));d=d0;
%   d.X=d.X(d.Y==1,:); d.Y=d.Y(d.Y==1,:); % only have one class available
%   [r a]=train(param(one_class_svm({kernel('rbf',.2),'optimizer="quadprog"'}),'nu', linspace(0.1,0.9,4)),d); 
%   subplot(411); plot(a{1},d);
%   subplot(412); plot(a{3},d);
%   subplot(413); plot(a{3},d);
%   subplot(414); plot(a{4},d);
%==============================================================================
% Reference : Estimating the support of a high-dimensional distribution
% Author    : B. Schölkopf , A. J. Smola , J. Platt , J. Shawe-Taylor and R. C. Williamson
% Link      : http://www.kernel-machines.org/papers/oneclass-tr.ps.gz
%==============================================================================
  
  %hyperparams 
  a.child=kernel;
  a.optimizer='libsvm';
  a.nu = 0;
  
  a.C = Inf;
 
  % AVAILABLE OPTIMIZERS (eventually):
  % andre,quadprog,libsvm 
      
  % model 
  a.alpha=[];
  a.alpha_cutoff=-1;
  a.b0=0;
  a.Xsv=[];
  
  p=algorithm('one_class_svm');
  a= class(a,'one_class_svm',p);
  a.algorithm.use_signed_output=1;
  
  if nargin==1,
    eval_hyper;
  end; 
  
 
