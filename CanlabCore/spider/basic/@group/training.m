function [r,a] =  training(a,d)
   
  if isa(d,'group') & strcmp(d.group,'one_for_each') & length(d.child)==length(a.child)

        %% pass only nth data object into nth algorithm

    r=[]; 
    for i=1:length(a.child)
      [res,a.child{i}]=train(a.child{i},d.child{i}); r{i}=res;
    end
	
  else  %% standard method
    
    r=[]; 
    for i=1:length(a.child)
      [res,a.child{i}]=train(a.child{i},d); r{i}=res;
    end
  end
  
  if length(r)==1 r=r{1}; else r=group(r); end;
  if isa(r,'group') r.group=a.group; end;
  



