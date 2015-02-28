
function a = chain(in,hyper)
 
%==================================================================
% chain object
%==================================================================
% A=CHAIN(I,H) returns a chain object initialized with 
%  a cell array of algorithms I and (hyper)parameters H.
%   
% This is used to create a chain of algorithms, the output of one
%  is fed into the input of the next.
%
% Examples: f=chain({fisher('output_rank=1;feat=5') knn});
%           [r,a]=train(cv({f knn}),toy);
%           get_mean(r)
%==================================================================
% Reference : 
% Author    : 
% Link      : 
%==================================================================
 
  if nargin == 0, in = {}; end
   
  if length(in)==1 a.child{1}=in; else a.child=in; end;
  
  % convert {} to group({})
  for i=1:length(a.child) 
     if isa(a.child{i},'cell') a.child{i}=group(a.child{i}); end;
  end
  
  p=algorithm('chain');
  a= class(a,'chain',p);
  %a.orig=a.child;  %% original algorithm set given (for retraining)
  
  %hyperparams
  if nargin==2
    eval_hyper;
  end
