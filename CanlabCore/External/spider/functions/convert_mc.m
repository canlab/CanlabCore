function [res]=convert_mc(ys,k)
  
% function [res]=convert_mc(ys,[k])
% 
%  convert_mc  - converts multi-class labels from examples x classes matrix
%                of {+1,-1}s to examples x 1 matrix of integers, and back
%                can provide optional k to say how many classes
%                there should be (incase some are missing)
%
%  e.g:
%  convert_mc([1 2 2]') returns [1 -1; -1 1 ; -1 1] and vice-versa
%  if you pass a data object instead, it converts the Y component
%  in the same way and returns the data object with converted Ys.
  
  if nargin==1 k=0; end;
  
  d=[];
  if isa(ys,'data')
    d=ys;
    ys=d.Y;
  end
  
  expand=1;
  if size(ys,2)==1 & min(min(ys))>0
    res=-ones(size(ys,1),max(max(ys),k));
    for i=1:size(ys,1) res(i,ys(i))=1; end;
  else   
    expand=0;
    if size(ys,2)==1
      res=ys; res(res==-1)=2;
    else
      for i=1:size(ys,1) res(i)=find(ys(i,:)==1); end;
      res=res';
    end
  end
  
  if ~isempty(d)
    d.Y=res;
    res=d;
  end
  


