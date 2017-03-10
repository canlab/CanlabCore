function dat = rmse(algo,dat)
  
  [x y]=get_xy(dat);
    
  if size(y,2)>1,  
      lss = 0;
      for i = 1:size(y,1)
        lss=lss + norm(x(i,:)-y(i,:))^2;
      end
      lss = lss/size(y,1);
  else
      lss = 0;
      for i = 1:size(y,1)
        lss=lss + norm(x(i)-y(i))^2;
      end
      lss = lss/length(y);
  end;
  lss = sqrt(lss);
  
  dat=data([get_name(dat) ' -> root mean squared error (2-norm)=' num2str(lss,4) ],[],lss);
