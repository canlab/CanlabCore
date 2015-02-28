
function [numOfEx,vecDim,outDim,numOfCls] = get_dim(dat)
  
%   [numEx,vDim,oDim,numCls] = get_dim(data)
%
%   Returns (ind given order) the number of examples, the vectors
%   dimension, the number of output dimensions and the number
%   of classes.
%
%   Note:
%   Sometimes the number of classes differs from the output dimension
%   (e.g. in binary pattern recognition often oDim equals 1 and numCls
%   equals 2).
  
  global X;
  global Y;
  
  numOfEx=length(dat.index); 
  vecDim=length(dat.findex);
  if isempty(dat.myY)
    outp=Y(dat.index,:);
    outDim=size(outp,2);
    numOfCls=outDim; 
    if numOfCls==1 
        fin1=find(outp==1); 
        fin2=find(outp==-1);
    end;
  else
    outp=dat.myY;
    outDim=size(outp,2); 
    numOfCls=outDim; 
    if numOfCls==1 
        fin1=find(outp==1); 
        fin2=find(outp==-1);
    end;
  end

  if numOfCls==1 
    if length(fin1)>0 & length(fin2)>0 & (length(fin1)+length(fin2))==numOfEx
      numOfCls=2;  %% -1s and +1s => classes
    end
  end
      




