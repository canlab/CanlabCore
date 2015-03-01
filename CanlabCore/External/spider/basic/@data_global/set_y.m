function [retDat] = set_y(dat,outp,ind,fInd)
 
%   data = set_y(data,Y)          sets the Y matrix of the data_global object
%   data = set_y(data,Y,indexes)  sets the Y matrix of the data_global object
%                                 for specific indexes
  
if nargin<3 | isempty(ind),
    ind=[1:size(outp,1)];
end
if nargin<4|isempty(fInd),
    fInd=[1:size(outp,2)];
end
  
if isempty(dat.myY)&nargin>2
    global Y;
    dat.myY=Y;
    disp('warning -- instantiated some indexes of y when using globals');
    keyboard
end
  
dat.myY(ind,fInd)=outp;
retDat=dat;
