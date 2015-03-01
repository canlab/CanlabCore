function [d,a] =  training(a,d)
  
  if a.algorithm.verbosity>0
    disp(['training ' get_name(a) '.... '])
  end
  
  [l n k k2]=get_dim(d); 
  d2=convert_mc(d,k2);
  d2=convert_mc(d2,k2);  %% try to get two columns even for binary problems
    
  if iscell(a.child)    %% make k copies of child alg, one for each class
    % do nothing
  else
    r=[]; for i=1:k r{i}=a.child; end;
    a.child=r; 
    a.prior=ones(1,length(a.child))/length(a.child);
  end
  
  for i=1:k 
    f=find(d2.Y(:,i)==1); dtmp=get(d2,f); 
    if(a.train_priors) a.prior(i)=length(f)/l; end;
    [r a.child{i}]=train(a.child{i},dtmp);
  end
  
  if ~a.algorithm.do_not_evaluate_training_error
    d=test(a,d);
  end
  
  
