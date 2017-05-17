function dat = roc50(algo,dat)
  
  
%
%           dat = roc50(algo,dat)
%
% calculates roc 50 curves. If param is not empty, a curve is also be
% plotted, where param determines the linestyle (e.g. param='--*')
%


  cutoff=50;

  [x y]=get_xy(dat);
    
  if sum(y==-1)==0  %% <-- not possible to produce a roc curve
    lss=0;
    r=data([get_name(dat) ' -> roc50 = ' num2str(lss,4) ],[],1);
    return;
  end

  [val ind]=sort(-x);
  x=-val; 
  y=y(ind);
  score=(y==1);
 
  numTot=length(score);
  numPos=sum(score==1);
  numNeg=sum(score==0);
  area=0; 
  height=0;  
  fps=0;

  for i=1:numTot
   if score(i)==1 
       height=height+1; 
   end;
   if score(i)==0 
       area=area+height; 
       fps=fps+1; 
       if cutoff<=fps 
           break; 
       end;
    end;
  end 
 
  lss=area/(fps*numPos);

  dat=data([get_name(dat) ' -> roc=' num2str(lss,4) ],[],lss);
  
  
