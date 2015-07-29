function dat = sensitivity(algo,dat)
  
  [x y]=get_xy(dat);
  if length(find(abs(x)==1))~=prod(size(x))
    disp('[warning: sign not taken before class_loss]');   
    x=sign(x); 
  end; 
  
  tPos=sum(y==1 & x==1);
  tNeg=sum(y==-1 & x==-1);
  fPos=sum(y==-1 & x==1);
  fNeg=sum(y==1 & x==-1);

  lss=tPos/(tPos+fNeg);
      
  dat=data([get_name(dat) ' -> sensitivity=' num2str(lss,4) ],[],lss);
