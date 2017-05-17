 function a = loss(comp,par,spr) 

%============================================================================ 
% Loss object - for calculating loss functions 
%============================================================================
% a=loss(lossType,param) 
%
%   calculates the difference between the input X and the ouput Y depending
%   on the specified loss type (for further information see below).
%   The loss can be calculated in two ways. The first is to it by training,
%   the second is to call the function with a data object as first
%   parameter (loss(d,loss_type,param)). The results are stored in the Y
%   part of the data object.
% 
%   Attributes (with defaults): 
%       type='class_loss'  -- type of loss (class_loss,linear_loss...)
%       param=[]           -- used parameters (can also be empty)
%
%   Methods:
%       train,test,calc
%
%   LOSS              |   PARAMETERS & DESCRIPTION
%   -------------------------------------------------------------------
%   class_loss        -- zero/one loss, L(x,y)=1 if x=y, 0 otherwise
%   confusion_matrix  -- matrix of [true-pos, false-pos; false-neg, true-neg]
%   epsilon_loss      -- L(x,y)= |x-y|, if |x-y|>epsilon, 0 otherwise
%   linear_loss       -- 1-norm, L(x,y)=|x-y|
%   one_class_loss    -- for one-class, e.g novelty detection, etc.
%   quadratic_loss    -- 2-norm, L(x,y)=|x-y|_2^2
%   roc               -- receiver/operator characteristic
%   roc50             -- receiver/operator characteristic, first n fps
%   sensitivity       -- tp/(tp+fn)
%   specificity       -- tn/(fp+tn)
%   kernel            -- loss derived from kernel matrix (param) 
%                     -- inner products in 'loss' space between examples.
%   alignment         -- L(x,y)= sum(sum( (x*x') .* (y*y'))) / normalization
%===================================================================================
% Reference : 
% Author    : 
% Link      : 
%===================================================================================
   
  %% <<--calculate loss like function (no objects is created)------>>
  if nargin>0 & isa(comp,'algorithm') 
    a.type='class_loss'; 
    if nargin>1 
        a.type=par; 
    end;
    
    a.param=[]; 
    if nargin==3, 
        a.param=spr; 
    end;
    
    algoTemp=algorithm(a.type);
    a= class(a,'loss',algoTemp);
    a=train(a,comp); 
    return;
  end 
  %% --------------------------------------------------------------------
   
   if nargin==0 
     a.type='class_loss';
   else 
     a.type=comp;
   end
   
   a.param=[];

   if nargin==2
     a.param=par;
   end;
   
   algoTemp=algorithm(a.type);
   a= class(a,'loss',algoTemp);
   
