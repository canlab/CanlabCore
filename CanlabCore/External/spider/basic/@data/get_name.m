function [nam] = get_name(dat)
 
%   get_name(DATA)      returns name describing the object
  nam=dat.algorithm.name;

%%%% SGE support
if isdeferred(dat), nam = get_name(dat.algorithm.deferred, nam, 1); end
%%%%
