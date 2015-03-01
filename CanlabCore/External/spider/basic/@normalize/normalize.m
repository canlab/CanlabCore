function a = normalize(hyper) 

%==============================================================================   
% Normalization object - simple normalization of data
%==============================================================================  
% A=NORMALIZE(H) returns a normalize object initialized with hyperparameters H. 
%
%  Performs different kinds of pre-processing on data.  
%
%  Hyperparameters, and their defaults
%  a.scale_type=0    -- 0: do nothing,   
%                       1: columns have mean 0, std 1
%                       2: rows have mean, std 1
%  		                3: try to do do both 1 & 2
%                       4: scale by correlation coefficients
%  a.sigmoid=0       -- yes or no, scale features nonlinearly through a
%                       sigmoid, i.e f(x)=tanh(sig*x) 
%
% 
%  Model parameters
%   mean_vec         -- threshold vectors  
%   scale_vec        -- scaling factors
%   sig              -- sigmoid steepness parameter
%   corr             -- correlation coefficients
%==============================================================================
% Reference : 
% Author    : 
% Link      : 
%==============================================================================

  %hyperparams
  a.scale_type=1;     % - 0: do nothing,   
                      %   1: columns have mean 0, std 1
                      %   2: rows have mean, std 1
		              %   3: try to do do both 1 & 2
                      %   4: scale by correlation coefficients
  a.sigmoid=0;        % - scale features nonlinear through a
                      %   sigmoid, i.e f(x)=tanh(s(x)) 
  
  % model 
  a.mean_vec=[];
  a.scale_vec=[];
  a.mean_vec2=[];
  a.scale_vec2=[];
  a.sig=[];
  a.corr=[];
   
  p=algorithm('normalize');
  a= class(a,'normalize',p);
 
  if nargin==1,
    eval_hyper;
  end;
