
function [res] = make_cell(dat)
% [res] = make_call(dat)   Converts a single data item into list of data items
  
  if ~isa(dat,'cell')
    res{1}=dat;
  else
    res=dat;
  end 
  