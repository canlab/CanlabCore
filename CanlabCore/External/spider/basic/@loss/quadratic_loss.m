function dat = quadratic_loss(algo,dat)
  
  [x y]=get_xy(dat);
  
  % Some regression algorithms give more outputs back than just the predictions.
  ydim=size(y,2);
  x=x(:,1:ydim);
    
  if size(y,2)>1,  
    lss=sum(sum(abs(x-y)').^2)/size(y,1);
  else
      lss=sum((x-y).^2)/length(y);
  end;
  
  dat=data([get_name(dat) ' -> quadratic_loss=' num2str(lss,4) ],[],lss);
