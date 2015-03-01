function [r,a] =  testing(a,d)
  
  if isa(d,'group') & strcmp(d.group,'one_for_each')
    
    %% PASS ONLY Nth DATA OBJECT INTO Nth ALGORITHM 
    
    r=[];
    for i=1:length(a.child)  
      r{i}=test(a.child{i},d.child{i});
    end
       
  else  %% STANDARD METHOD 
    r=[];
    for i=1:length(a.child)
      r{i}=test(a.child{i},d);
    end
  end
    
  if length(r)==1 r=r{1}; else r=group(r); end;
  if isa(r,'group') r.group=a.group; end;
