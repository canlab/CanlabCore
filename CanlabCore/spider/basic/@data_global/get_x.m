function [retX] = get_x(dat,ind,fInd)
%   get_x(data)         returns X matrix of the data_global object
%   get_x(data,indexes) returns X matrix of the data_global object 
%                        for given indexes
global X;
   
if nargin<2|isempty(ind),
  if isempty(dat.myX)   
      ind=[1:length(dat.index)];  
  else                
      ind=[1:size(dat.myX,1)];    
  end;
end;  
if nargin < 3 | isempty(fInd),
  if isempty(dat.myX)   
      fInd=[1:length(dat.findex)];  
  else                
      fInd=[1:size(dat.myX,2)];    
  end;
end
  
if ~isempty(dat.myX)  
  retX=dat.myX(ind,:);
else
  retX=X(dat.index(ind),dat.findex(fInd)); 
end
 
