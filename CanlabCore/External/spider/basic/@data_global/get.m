function [dat] = get(d,indx,findx)
 
%   data=get(data)                  returns the data of the data_global object  
%   data=get(data,index,featInd)    returns the data for given indexes for
%                                   examples and features 
 
empty=0; if nargin>1 empty=isempty(indx);  end;
if nargin==1 | empty
  if isempty(d.myY)   indx=[1:length(d.index)];  
  else                indx=[1:size(d.myY,1)];    end;
end;  
        
if nargin<3 
  if isempty(d.myX)   findx=[1:length(d.findex)];  
  else                findx=[1:size(d.myX,2)];    end;
end;  
 
xx = d.myX;
if ~isempty(xx)
  xx=xx(indx,findx);
end
yy = d.myY;
if ~isempty(yy)
  yy=yy(indx,:);
end
inds  = d.index(indx);
finds = d.findex(findx);
dat   = data_global(d.algorithm.name,xx,yy,inds,finds); 
 
