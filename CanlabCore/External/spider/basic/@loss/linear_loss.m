function dat = linear_loss(algo,dat)
  

  [x y]=get_xy(dat);
  
  if size(y,2)>1,  
    loss=sum(sum(abs(x-y)'))/size(y,1);
  else
      loss=sum(abs(x-y))/length(y);
  end;
  
  dat=data([get_name(dat) ' -> linear_loss='  num2str(loss,4) ],[],loss);
