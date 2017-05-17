function [res,a] =  training(a,dat)
  
% [results,algorithm] =  train(algorithm,data,trn)
  
  disp(['training ' get_name(a) '.... '])
  
  [numEx,vDim,oDim] = get_dim(dat);
  a.nrofclasses=oDim;
  yTemp=get_y(dat);
 
  lent = length(a.child); 
  childIndex=0;
  for i=1:oDim,
    for j=(i+1):oDim,
      childIndex=childIndex+1;
      if (childIndex>lent),
        a.child{childIndex} = a.child{mod(childIndex,lent)+1};
      end;
      indTemp_i = find(yTemp(:,i)==1);
      indTemp_j = find(yTemp(:,j)==1);
      indTemp = [indTemp_i',indTemp_j'];
      leni = length(indTemp_i);
      lenj = length(indTemp_j);
      Yij = ones(leni+lenj,1);
      Yij(leni+1:lenj+leni)=-1;
      datTemp = get( dat, indTemp);
      datTemp = set_y( datTemp, Yij);
      %  datTemp = data('tmp',xTemp,Yij);    %%%% looses indices of data!!!
      datTemp=set_name(datTemp,['Machine ' num2str(childIndex)]);
      [r{childIndex},a.child{childIndex}]=train(a.child{childIndex},datTemp);
    end;  
  end;
  
  res=test(a,dat);
