
function a = rss_mp(alg,hyper) 

%=================================================================
% RSC_MP reduced set selection
%=================================================================
% a = rss_mp(alg,hyper)
% generates a rss object, using the matching pursuit selection method
%
% hyperparameters:
% child=svm         algorithm worked on
% max_loops=1e5     maximum number of basis functions
% tolerance=1e-5    tolerated loss in ||w-w*||^2
% backfit=1         backfit on every nth iteration
% backfit_at_start=100  always backfit for first e.g. 100 iterations 
% dont_revisit=1    dont return to old basis function optimization, always get new one
% reoptimize_b=1    recalculate the threshold b0
% alpha_cutoff=0    throw away svs with abs(alpha)<n
% bal_w=0           treat multiple w's as equal by normalizing by length
% optimizer='iterative' iterative update of matrix inverse
% 
% model:
% alpha         new alphas for rs-vectors
% Xsv           rs vectors
% b0            the threshold
%
% stats:
% w2=0          final value of ||w-w*||^2 
% res=[]        results on a separate test set
% dtst=[]       separate test set  
% test_on=0     iterations to test on
%
% methods:
% train         constructs a reduced set, returns trained rs-machine
% test          tests new rs-machine on supplied data
%
% example:
% d=gen(toy2d('2circles','l=100'));
% [r,a]=train(svm({kernel('rbf',1),'C=10000','alpha_cutoff=1e-2'}),d);
% [r,a2]=train(rss_mp(a,'tolerance=1e-2'),d);
% test(a2,d,loss)
%
%=================================================================
% author: goekhan bakir, jason weston
% reference: fast binary and multi-output rss, 2004
%=================================================================

  
  %hyperparams 
  a.child=svm;
  a.max_loops=100000;% maximum number of basis functions
  a.tolerance=1e-5; 
  a.backfit=1;      % backfit on every nth iteration
  a.backfit_at_start=100; % always backfit for first e.g. 100 iterations 
  a.dont_revisit=1;  % dont return to old basis function optimization, always get new one
  a.reoptimize_b=1;
  a.alpha_cutoff=0;
  a.bal_w=0;         % treat multiple w's as equal by normalizing by length
  a.optimizer='iterative'; % iterative update of matrix inverse
  
  % model 
  a.alpha=[];
  a.Xsv=[];  
  a.b0=0;
  
  % training / testing statistics
  a.w2=0;  % final value of ||w*-w^2||^2
  a.res=[];% results on a separate test set
  a.dtst=[];% separate test set  
  a.test_on=0; % iterations to test on
  
  if nargin==0
    a.child=svm;  
  else 
    a.child=alg; %% algorithm to use  
  end
  
  p=algorithm('rss_mp');
  a= class(a,'rss_mp',p);
  a.algorithm.use_signed_output=0;
  
  if nargin==2
    eval_hyper;
  end  
  
 
