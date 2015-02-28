function [retX,retY] = get_xy(dat,ind)
 
%   get_xy(data)           returns X and Y matrices of the data_global object
%   get_xy(data,indexes)   returns X and Y matrices of the data_global object for given
%                           indexes
  
    
global X; 
global Y;

if nargin==1 
  if isempty(dat.myX)   
      ind=[1:length(dat.index)];  
  else                
      ind=[1:size(dat.myX,1)];    
  end;
end;  
if ~isempty(dat.myX)  
  retX=dat.myX(ind,:);
else
  retX=X(dat.index(ind),dat.findex); 
end
if ~isempty(dat.myY)  
  retY=dat.myY(ind,:);
else
  retY=Y(dat.index(ind),:); 
end
