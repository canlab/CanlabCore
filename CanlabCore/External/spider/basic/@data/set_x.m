function [retDat] = set_x(dat,inp,ind,fInd)
 
%   DATA=set_x(DATA,X)         sets the X matrix of data-object
%   DATA=set_x(DATA,X,INDEXES) sets the X matrix of data-object
%                              for given indexes
    
if nargin==2 | (isempty(ind)&isempty(fInd)),
  dat.X=inp;
  if ~(size(inp,1)==length(dat.index))
    dat.index  = [1:size(inp,1)];
  end
  if ~(size(inp,2)==length(dat.findex) )
    dat.findex = [1:size(inp,2)];
  end
else
    if isempty(ind),
        ind = [1:size(inp,1)];
    end;
    if (nargin<4) | isempty(fInd),
        fInd = [1:size(inp,2)];
    end;
    dat.X(ind,fInd)=inp;
end
retDat=dat;