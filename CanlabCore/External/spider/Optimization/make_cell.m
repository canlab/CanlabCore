
function [res] = make_cell(d)
  %% converts single data item into list of data items
  
  if ~isa(d,'cell')
    res{1}=d;
  else
    res=d;
  end 
  