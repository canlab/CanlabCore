function [retVal]=get_norm(k,dat,ind)
% [val]=get_norm(kernel,data,[index])
%  Returns the norm of the vectors via the data and
%   kernel specified, optionally one can calculate this only
%   for some of the examples using index
  
if k.calc_on_output %% output kernels
    dat.X=dat.Y; 
end  

if nargin==2
    tmp = get_dim(dat);
    ind=[1:tmp(1)];
end
if strcmp(k.ker,'linear'),
  tmp = get_x(dat,ind);  
  retVal = sqrt(sum(tmp.^2,2));
else
% slow implementation but it works
  %for i=1:length(ind),
  %	  retVal(i) = get_kernel(k,get(dat,ind(i)));
  %end;
  retVal = get_kernel(k,get(dat,ind)); 
  retVal = diag(retVal)';
  retVal = sqrt(retVal');
end;



