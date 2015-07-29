function dat = balanced_loss(algo,dat)
  
  [x y]=get_xy(dat);
    
  if isempty(find(abs(x)==1)) 
    disp('[warning: sign not taken before balanced_loss]');
    x=sign(x); 
  end; 
    
  if size(y,2)>1,  
      error('not implemented for multi-class yet')
  else
    if sum(y==1)>0
      lss1= sum(x==1  & y==1) / sum(y==1);
    else lss1=1; end;
    if sum(y==-1)>0
      lss2= sum(x==-1 & y==-1) / sum(y==-1);
    else lss2=1; end;
    lss=1-(lss1+lss2)/2;
  end;
  
  dat=data([get_name(dat) ' -> balanced_loss=' num2str(lss,4)],[],lss);
