function [nam] = get_name(dat)
 
%   get_name(data)  returns the name of the object
  nam=dat.algorithm.name;

%%%% SGE support
if isdeferred(dat), nam = get_name(dat.algorithm.deferred, nam, 1); end
%%%%
