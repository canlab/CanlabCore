function [d,a] =  testing(a,d)
        
  for i=1:length(a.child)  
    d=test(a.child{i},d); 
  end
    