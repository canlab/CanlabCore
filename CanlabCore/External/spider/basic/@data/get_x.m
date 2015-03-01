function [retX] = get_x(dat,ind,fInd)
%   get_x(DATA)          returns X matrix of the data object
%   get_x(DATA,INDEXES)  returns X matrix of the data object 
%                        for given indexes
  
if (nargin==1)|isempty(ind),
    ind=[1:size(dat.X,1)]; % <---- return all we got 
end;    
if nargin < 3 | isempty(fInd),
    fInd = [1:size(dat.X,2)];
end;
if ~isempty(dat.X)  
  retX=dat.X(ind,fInd);
else
  retX=[]; 
end
