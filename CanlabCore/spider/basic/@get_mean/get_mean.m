 function [a,s] = get_mean(d,l,hyper) 

%==========================================================================
% get_mean object - creates object to calculate mean of groups
%==========================================================================
%  R=get_mean(A,[L])
% 
%  An object which calculates the mean(s) R of results with algorithms A,
%  together with the standard error (std. dev / sqrt(number of trials))
%  optionally calculating the loss with loss type L (see help loss)
%  e.g 'class_loss',.. (NOTE: input is a string, not a loss object).
%  when trained, e.g:  train(get_mean(cv(svm)),gen(toy))
%
%  Hyperparameters:
%   take_average=0       -- method of taking average over groups
%                           1: average over inside leaf of tree of groups
%                           2: average over outside leaf of tree of groups
%                           0: try to guess which is most appropriate 
%
%        'take_average' is necessary because if you have a group
%         of groups of algorithms you could wish to take the mean in two
%         different ways: e.g with group({group({a b} group({c d}))}) 
%         you may either wish to average a&b and c&d or a&c and b&d.
%
%         NOTE: The same effect can usually be achieved by doing
%               get_mean(r) and get_mean(r') with a data object r
%   
%  Alternatively, can be used in a feed
%  forward network, e.g: train(chain({ group(cv(svm)) get_mean }),gen(toy))
%  
%  Can also be called like a function with 
%  A=get_mean(D,[L]) calculates the mean for groups of data objects D.
%=========================================================================
% Reference : 
% Author    : 
% Link      : 
%=========================================================================
  
  %% --------- don't make object : calculate get_mean like function --------
  if nargin>0 & isa(d,'algorithm') %% input data - return calculated get_mean
    if d.is_data

      if nargin<3 hyper=[]; end; a=[]; eval_hyper;
      take_average=0; 
      if isfield(a,'take_average') take_average=a.take_average; end;
 
      if nargin==1 [a s]=calc_mean(d,[],0,take_average); 
         else [a s]=calc_mean(d,l,0,take_average); end;
      return;
    end
  end 
  %% --------------------------------------------------------------------
  
  a.loss_type=[]; a.child=[];
  a.take_average='guess';
  if nargin>0    
    if isa(d,'algorithm') 
      a.child=d;
      if nargin>1 a.loss_type=l; end;
    else
      if nargin>0 a.loss_type=d; end;
    end
  end;
  p=algorithm('get_mean');
  a= class(a,'get_mean',p);
   
  if nargin<3 hyper=[]; end;
  eval_hyper;
