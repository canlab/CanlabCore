function [retY] = get_y(dat,ind,fInd)
 
%   get_y(data)     returns Y matrix of the data_global object
%   get_y(data,indexes) returns Y matrix of the data_global object for given
%                   indexes
  
global Y;
     
if nargin<2|isempty(ind),
  if isempty(dat.myY)   
      ind=[1:length(dat.index)];  
  else                
      ind=[1:size(dat.myY,1)];    
  end;
end;  
if nargin < 3 | isempty(fInd),
  if isempty(dat.myY)   
      fInd=[1:size(Y,2)];  
  else                
      fInd=[1:size(dat.myY,2)];    
  end;
end
  
if ~isempty(dat.myY)  
  retY=dat.myY(ind,fInd);
else
  retY=Y(dat.index(ind),fInd); 
end
