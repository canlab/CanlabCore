function [dat] = get(dt,ind,fInd) 

%   data=get(data)                      returns the data of the data object 
%   data=get(data,examIndex,featIndex)  returns the data of the data object 
%                                       for given indexes   

empty=0; 
if nargin>1 
    empty=isempty(ind);  
end;

if (nargin==1 | empty) 
    ind=[1:size(dt.Y,1)];  %<--- return all examples 
end;  
if nargin<3 
    fInd=[1:size(dt.X,2)]; %<--- return all features
end;          

y=[]; 
if ~isempty(dt.Y) 
    y=dt.Y(ind,:); 
end;

x=[]; 
if ~isempty(dt.X) 
        x=dt.X(ind,fInd); 
end;
  
dat  = data(dt.algorithm.name,x,y); 
dat.index = dt.index(ind);
dat.findex = dt.findex(fInd);



