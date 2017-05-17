
function [r2,mn,st]= cv_mean(a,r,loss_type,param1)
  
if nargin<3
  [r2,st]=get_mean2(r);
  [mn,st]=get_mean2(r,[],1);
else
  [r2,st]=get_mean2(r,loss_type);
  [mn,st]=get_mean2(r,loss_type,1);
end