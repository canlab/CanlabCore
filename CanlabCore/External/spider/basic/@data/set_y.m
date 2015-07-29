function [retDat] = set_y(dat,outp,ind,fInd)
 
%   data=set_y(DATA,Y) sets the Y matrix of data
%   data=set_y(DATA,Y,INDEXES) sets the Y matrix of data for given indexes
  
if nargin==2 | (isempty(ind)&isempty(fInd)),
  dat.Y=outp;
else
    if isempty(ind),
        ind = [1:size(outp,1)];
    end;
    if (nargin<4)|isempty(fInd),
        fInd = [1:size(outp,2)];
    end;
    dat.Y(ind,fInd)=outp;
end
retDat=dat;
