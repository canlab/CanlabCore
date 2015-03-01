function [retK]=  dist(algo,x,y)
num=size(x,1);
outDim=size(y,2);
retK=diag(x*x')*ones(1,outDim) +  ones(num,1)*diag(y'*y)' -2*(x*y);
