function [l,n,k,k2] = get_dim(dat)  
%   [numEx,vDim,oDim,numCls]=get_dim(data)
%    return (in given order) the number of examples, the dimension of vectors
%    the number of output dimensions and the number of classes.
%    
%    Note:
%    Sometimes the number of classes differs from the output dimension
%    (e.g. in binary pattern recognition often oDim equals 1 and numCls
%    equals 2).



  l=length(dat.index);
  n=length(dat.findex);
  k=size(dat.Y,2); 
  k2=k;
  
  if k2==1
    f1=find(dat.Y==1); f2=find(dat.Y==-1); 
    if length(f1)>0 & length(f2)>0 & (length(f1)+length(f2))==l
      k2=2;  % -1s and +1s => two classes
    end
  end