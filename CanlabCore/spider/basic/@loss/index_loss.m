function dat = index_loss(algo,dat)
  
  [x y]=get_xy(dat);
  res = sum((x-y).^2,2);
  lss = length(find(res > 0))/max(size(res));  
  
  dat=data([get_name(dat) ' -> index_loss=' num2str(lss,4) ],[],lss);
