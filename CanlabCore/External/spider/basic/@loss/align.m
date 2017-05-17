function dat = align(algo,dat)
 
  [x y]=get_xy(dat);
  
  if size(x,1)~=size(x,2) 
      x=x*x'; 
  end;
  if size(y,1)~=size(y,2) 
      y=y*y'; 
  end;
    
  loss=sum(sum( x .* y)) / sqrt(sum(sum(x.^2)) * sum(sum(y.^2)) );
   
  dat=data([get_name(dat) ' -> alignment='  num2str(loss,4) ],[],loss);

