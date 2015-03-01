function ret = epsilon_loss(algo,dat,ind)
  
  if nargin==2          %% no indexing, use all data
    [x y]=get_xy(dat);
  else
    [x y]=get_xy(dat,ind);
  end
  
  tmp = abs(x-y) > algo.param;
  
  if size(y,2)>1,  
      loss=sum(sum(tmp.*abs(x-y)))/(size(y,1)); 
  else
      loss=sum(tmp.*abs(x-y))/length(y);
  end;
  
  ret=data([get_name(dat) ' -> epsilon_loss='  num2str(loss,4) ],[],loss);