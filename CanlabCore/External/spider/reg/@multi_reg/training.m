function [res,a] =  training(a,d)
  
% [results,algorithm] =  train(algorithm,data,trn)
  

  disp(['training ' get_name(a) '.... '])

 [m,n,Q] = get_dim(d);
 a.Q=Q; % store number of outputs for testing
 dtmp=d; 
 Ytmp=get_y(dtmp);
 
  for i=1:Q,
    if (i>length(a.child)),
      a.child{i} = a.child{mod(i,length(a.child))+1};
    end;
    dtmp=set_y(dtmp,Ytmp(:,i));
    dtmp=set_name(dtmp,['Machine ' num2str(i)]);
    [r{i},a.child{i}]=train(a.child{i},dtmp);
  end;
  
  res=test(a,d);


  

