function [a]=add(a,d)
  
% function [a]=add(a,d)
%
% add a single algorithm or a cell array algorithms to an algs
% object
  
  d=make_cell(d);
  r=a.child;
  if ~isempty(r) r={r{1:length(r)} d{1:length(d)}}; else r=d; end
  a.child=r;