function a = algorithm(c) 
  

%=============================================================================================
%   ALGORITHM algorithm object      
% ============================================================================================
%           a=algorithm(hyperParam) 
%
%   An algorithm with given hyperparamters is constructed.
%   All algorithms will inherit these hyperparamters resp. option settings.
%
%
%   Hyperparameters/Option Settings (with defaults):
%       training_time=0                   -- cputime needed for training is taken
%       trained=0                         -- true, if a model has been learnt
%       do_not_retrain=0                  -- determines if the algorithm has to
%                                            be retrained (useful for param)
%       use_prev_train=0                  -- determines if previous model shall
%                                            be used (e.g. param objects don't need retraining)
%       do_not_evaluate_training_error=0  -- selfexplanatory (speed up computation)
%       use_signed_output=1               -- true if sign after output
%                                            shall be taken
%       verbosity=1                      -- verbosity level 
%       is_data=0                        -- true if object is considered as data 
%       alias=[]                         -- alternative names for this
%                                           object
%                                        -- accessing members, e.g can do:
%                                           a=svm;
%                                           a.alias={'p1','C','p2','alpha'}; 
%                                           a.p1=5; a.p2=ones(1,10);
%       deferred = []; 
%
% Methods:
%  train,test                            -- selfexplanatory
%=============================================================================================
% Reference : 
% Author    : 
% Link      : 
%=============================================================================================
  
  if nargin==0  
    a.name='algorithm';
  else
    a.name=[c]; 
  end

  %%<<-----calculate, if need to retrain when trying new parameters ----->
  a.trained=0;
  a.do_not_retrain=0;
  a.training_time=0;
   %%     ---used when training takes into account the current params
   %%        of the algorithm to train faster e.g.: with features, or testing
   %%        a new C,...
  a.use_prev_train=0;
  
  %%%%%%--do not calculate training error (for speeding up)
  a.do_not_evaluate_training_error=0;
  
  %%%%%%---wheter to use signed output in classifier or not
  a.use_signed_output=1; 
  
  %%%%%%---use to define verbosity of output
  a.verbosity=1; 
  
  %%%%%%---use to tell if data or not (e.g for objects like wilcox) 
  a.is_data=0;

  %%%%%%---alternative names for accessing members
  a.alias=[];
  
  a.deferred = [];
  
  a= class(a,'algorithm');
  
  




