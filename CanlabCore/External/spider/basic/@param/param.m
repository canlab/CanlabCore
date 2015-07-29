function a = param(alg,param,values,hyper) 
%===========================================================================         
% PARAM param object
%===========================================================================
% A=PARAM(I,N,P,H) returns a param object which
%  enumerates the possible values P of hyperparameter N of algorithm I
%  and with (hyper)parameters H.
%
% This is used to create a set of algorithms with different hyperparameters
%  e.g k-NN for different values of k.
%
% Example:  f=param(knn,'k',[1:2:10]); [r,a]=train(cv(f),generate(toy));
%           get_mean(r)
%
% Can also be used to enumerate multiple parameters by specifying N
% and P as cell arrays, e.g:
%  train(param(svm(kernel('poly',2)),{'C','kerparam'},{[1 100],[1 2 3]} ),toy)
% data=hyper(algorithm,parameter_name,parameter_values)
%===========================================================================
% Reference : 
% Author    : 
% Link      : 
%===========================================================================
  % runs a given algorithm many times with different values
  % of a given parameter    
% model
  if nargin == 0 % need to allow this possibility because matlab sometimes calls constructors silently without input args
  	alg = [];
	param = [];
	values = [];
  end
  
  a.child=alg;
  a.group='separate'; %% grouping of output data objects
  a.param=param;
  a.values=values;
  a.store_all=1;      %% store all objects after param search
  a.force_no_train=0; %%  can force to use model already learnt 
                      %% (e.g change value of b0, rfe feat, etc)
  a.force_use_prev_train=0; %% can force the models to use their previous
                            %% values so that it trains faster
			    
  par=algorithm('param');
  a= class(a,'param',par);
%  if nargin < 3, error('three input arguments required'), end
  if nargin==4,
    eval_hyper;
  end;  
  
  
  if a.store_all  %% ---- enumerate all possibilities
    a.child=[];
    if ~isa(a.values,'cell') val{1}=a.values; else val=a.values; end
    if ~isa(a.param,'cell') p{1}=a.param; else p=a.param; end
    for i=1:length(val) sz(i)=length(val{i}); end;
    tot=prod(sz); %% find permutations of hyperparameters
    
    for i=1:tot %% all permutations
      vars=num2choice(a,i,sz); %% hyperparams for this iteration      
      for j=1:length(p)      %% set hypers       
	v=val{j}(vars(j));  
	eval(['alg.' p{j} '=v;']); 
	
	%% set whether to force/ not force training for that object 
	e=p{j}; t=max(find(e=='.')); f=a.force_no_train;
	eval(['alg.' e(1:t) 'algorithm.do_not_retrain=f;'  ]);
      end
      a.child{i}=alg;
    end
  end
