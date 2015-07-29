function [dat] = set_x(dat,inp,ind,fInd)
 
%   data=set_x(data,x) sets the X matrix
%   data=set_x(data,x,ind) sets the X matrix for specific indexes

if nargin==2 | (isempty(ind)&isempty(fInd)),
  dat.myX=inp;
  if ~(size(inp,1)==length(dat.index))
    dat.index  = [1:size(inp,1)];
  end
  if ~(size(inp,2)==length(dat.findex) )
    dat.findex = [1:size(inp,2)];
  end
  return;
end
     
if nargin<3 | isempty(ind),
    ind=[1:size(inp,1)];
end
if nargin<4|isempty(fInd),
    fInd=[1:size(inp,2)];
end
  
if isempty(dat.myX)&nargin>2
    global X;
    dat.myX=X;
    disp('warning -- instantiated some values of x when using globals');
end
  
dat.myX(ind,fInd)=inp;








