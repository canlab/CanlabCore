function dat = N_class_loss(algo,dat)
  
  x = dat.X;
  y = dat.X;
%   [1:length(x)]';
  
  lss = x-y;
  lss(lss ~= 0) = 1;
  lss = sum(lss)/length(x);
  
  dat=data([get_name(dat) ' -> N_class_loss=' num2str(lss,4) ],[],lss);

