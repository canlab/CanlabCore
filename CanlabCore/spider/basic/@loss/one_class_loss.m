function dat = one_class_loss(algo,dat)
  
  [x y]=get_xy(dat);
  
  y = ones(size(x)); 
  lss=sum((x~=y))/length(y);
   
  dat=data([get_name(dat) ' -> one_class_loss='  num2str(lss,4) ],[],lss);
