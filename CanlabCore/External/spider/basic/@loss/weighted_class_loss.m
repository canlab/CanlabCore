function dat = weighted_class_loss(algo,dat)
  
  [x y]=get_xy(dat);
  if length(find(abs(x)==1))~=prod(size(x))
    disp('[warning: sign not taken before class_loss]');   
    x=sign(x); 
  end; 
  W=algo.param;
  if size(y,2)>1,  
    sz1=size(x,2); 
    sz2=size(y,2); 
    if sz1~=sz2 
      tmp=-ones(size(x,1),max(sz1,sz2)); 
      tmp(:,1:sz1)=x; 
      x=tmp;
      y=-ones(size(y,1),max(sz1,sz2)); 
      tmp(:,1:sz2)=y; 
      y=tmp; 
    end;
    lss=sum(W'*(x~=y))/(2*size(y,1)); 
  else
    lss=(W'*(x~=y))/length(y);
  end;
  
  dat=data([get_name(dat) ' -> class_loss=' num2str(lss,4) ],[],lss);

