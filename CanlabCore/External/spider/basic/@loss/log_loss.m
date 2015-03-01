function dat = log_loss(algo,dat)
  
  [x]=get_x(dat);
  
  loss= sum(-log(x))/size(x,1);
    
  dat=data([get_name(dat) ' -> log_loss=' num2str(loss,4) ],[],loss);

