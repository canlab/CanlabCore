function dat = BER(algo,dat)
  
  [x y]=get_xy(dat);     
  if length(find(abs(x)==1))~=prod(size(x))
    disp('[warning: sign not taken before class_loss]');   
    x=sign(x); 
  end; 
  
  if size(y,2)>1
      error('not implemented for multi-class yet')
  else
	tPos=sum(y==1 & x==1);  %%  d
	tNeg=sum(y==-1 & x==-1); %% a
	fPos=sum(y==-1 & x==1); %% b
	fNeg=sum(y==1 & x==-1); %% c

	lss=0.5*(fNeg/(tPos+fNeg) + fPos/(tNeg + fPos));
      
  	dat=data([get_name(dat) ' -> BER=' num2str(lss,4) ],[],lss);

  end;