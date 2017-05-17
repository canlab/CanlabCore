
function [d] = concatenate(tres,res,trn,tst)

%  function [d] = concatenate(trainRes,testRes,train,test)
%
%   Returns unified data (e.g with train and test results)
%   with trnRes in specified indexes train and tstRes in 
%   indexes test.
      
  global X; global Y;
 
  [l1,n1,k1]=get_dim(res);
  [l2,n2,k2]=get_dim(tres);
  l=l1+l2; 
    
  globalx=1; globaly=1;
  if ~isempty(tres.myX) | ~isempty(res.myX)  globalx=0; end
  if ~isempty(tres.myY) | ~isempty(res.myY)  globaly=0; end
    
  d=data_global(res.name);
  d.findex=res.findex; % <--- assumption: tres.findex = res.findex
  d.index=ones(1,l);   % <--- total indexes = tres+res
  if ~isempty(trn) d.index(trn)=tres.index; end;
  if ~isempty(tst) d.index(tst)=res.index; end;
  if ~globalx
    Xs=zeros(l,n2); 
    if size(Xs,2)>0
    if ~isempty(trn) Xs(trn,:)=get_x(tres);  end;
    if ~isempty(tst) Xs(tst,:)=get_x(res); end;
    end
    d=set_x(d,Xs);        %% <--- store results into a data object     
 end
  if ~globaly 
    Ys=zeros(l,size(tres.Y,2));   
    if ~isempty(trn) Ys(trn,:)=get_y(tres); end; 
    if ~isempty(tst) Ys(tst,:)=get_y(res); end;
    d=set_y(d,Ys);
  end
