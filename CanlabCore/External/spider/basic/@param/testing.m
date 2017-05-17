function [r] =  testing(a,d)
  
% [results,algorithm] =  train(algorithm,data)
  if a.store_all==0
    error('cannot run test on param object when store_all=0')
  end
  
  
  if ~isa(a.values,'cell') val{1}=a.values; else val=a.values; end
  if ~isa(a.param,'cell') p{1}=a.param; else p=a.param; end
   
  for i=1:length(val) sz(i)=length(val{i}); end;
  tot=prod(sz); %% find permutations of hyperparameters 
  
  r=[];  
  for i=1:tot %% all permutations    
    
    if isa(d,'group') & strcmp(d.group,'one_for_each') & length(d.child)==length(a.child)
      res=test(a.child{i},d.child{i});
    else
      res=test(a.child{i},d);
    end
    r{i}=res;
  end
  
  if length(r)==1 r=r{1}; else r=group(r); end;
  if isa(r,'group') r.group=a.group; end;
  

