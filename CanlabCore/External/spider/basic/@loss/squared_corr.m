function dat = squared_corr(algo,dat)
  
  [x y]=get_xy(dat);
 
  lss=corrcoef(x,y);
  lss = lss.*lss;
  l = lss(1,2);  
  dat=data([get_name(dat) ' -> squared_corr=' num2str(l,4) ],[],l);
