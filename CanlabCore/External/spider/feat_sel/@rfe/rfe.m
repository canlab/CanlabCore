function a = rfe(c,hyper)     

%=========================================================================
% RFE recursive feature elimination object
%=========================================================================   
% A=RFE(C,[H]) returns an rfe object initialized with base algorithm
% C (usually an svm) and hyperparameters H. 
%
% RFE selects features via greedy backward selection: it trains a
% base algorithm, and calls the method get_w of the base algorithm
% and removes the feature(s) with the smallest absolute value returned
% from this method. get_w in a linear svm for example returns the
% weights for each dimension defining the hyperplane. Base algorithms
% with a method get_w include: one_vs_one, one_vs_rest, lda, stab, svm, svr
% 
%  Hyperparameters, and their defaults
%
%   feat=[]        -- number of selected features to be selected, if feat=0, 
%                     then the minimum number of feature will be computed.
%   feature_group  -- Use this to supply a vector containing a group
%         = []        index for each feature. Features with the same
%                     group index are eliminated together. By default,
%                     the vector 1:n is used, meaning features are
%                     eliminated independently.
%   force=[]       -- one value per feature or feature group. Features/groups
%                     with negative values will be eliminated first (in
%                     increasing order of force value). Features/groups
%                     with positive values will be eliminated last (again,
%                     in increasing force value order). A value of 0 means
%                     that feature elimination order is determined by the
%                     normal method.
%   speed=0        -- rfe only removes a single feature each iteration 
%                     if less than (speed) features left, otherwise it
%                     halves the number of features each iteration.
%   output_rank=0  -- whether a ranking is desired, if set to 0 then a
%                     classification is perfomed after feature selection.
%   child=svm      -- The base classifier to be used at each step
%
%  Model
%
%  a.rank=[]       -- ranking of features
%  a.child=svm     -- base classifier trained at end of process
%
%  Methods:
%   train, test
%
% Example:
% d=gen(toy); a=rfe; a.feat=10; a.output_rank=1;[r,a]=train(a,d);
% a.rank  % - lists the chosen features in  order of importance
%
%=========================================================================
% Reference : Gene selection for cancer classification using support vector machines
% Author    : I. Guyon, J. Weston, S. Barnhill, and V. Vapnik
% Link      : http://www.stanford.edu/class/cs374/CLUSTER_Guyon.pdf
%========================================================================= 

% data=rfe(algorithm)
  % runs rfe on a given algorithm (default=linear binary svm)
  % requires an algorithm which has a function get_w
    
  %hyperparams
  a.feat=[];       % number of features (or feature groups, see feature_group)
  a.speed=0;       % rfe only removes a single feature if less than
                   % speed features left 
  a.output_rank=0; % output labels, not selected features 
    
  
  a.feature_group = []; % Use this to supply a vector containing a group
                        % index for each feature. Features with the same
                        % group index are eliminated together. By default,
                        % the vector 1:n is used, meaning features are
                        % eliminated independently.

  a.force = [];
  
  a.test_on_the_fly.loss =loss('class_loss');
  a.test_on_the_fly.data =[]; % cv/training sets this, allowing evaluation
                              % of cv error at each feature elimination---
							  % otherwise you would have to go back and
							  % re-train the classifier once for every
							  % feature count you want evaluated, which
							  % doubles the amount of time taken.
												
						
  % model
  a.rank=[];
  
  if nargin==0
    a.child=svm;  
  else 
    a.child=c; %% algorithm to use  
  end
  
  p=algorithm('rfe');
  a= class(a,'rfe',p);
  if nargin==2
    eval_hyper;
  end  
  
