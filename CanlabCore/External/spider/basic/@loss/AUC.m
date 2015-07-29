function dat = AUC(algo,dat)
  
  [x y]=get_xy(dat);     
  if length(find(abs(x)==1))~=prod(size(x))
    disp('[warning: sign not taken before class_loss]');   
    x=sign(x); 
  end; 
  
  if size(y,2)>1
      error('not implemented for multi-class yet')
  else
	posidx=find(y>0);
	negidx=find(y<0);
	[p1,p2]=size(posidx);
	[n1,n2]=size(negidx);
	posout=repmat(x(posidx),n2,n1);
	negout=repmat(x(negidx)',p1,p2);
	rocmat=2*(negout<posout);
	rocmat(negout==posout)=1;
	lss=sum(sum(rocmat))/(2*max(n1,n2)*max(p1,p2));	
      
  	dat=data([get_name(dat) ' -> AUC=' num2str(lss,4) ],[],lss);

  end;