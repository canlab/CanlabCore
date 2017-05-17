function [d,a] =  training(a,d)
        
  for i=1:length(a.child)  
    [d a.child{i}]=train(a.child{i},d);
    if i<length(a.child) & (isa(a.child{i},'group') | isa(a.child{i},'param')) 
      a.child{i}.group='one_for_each'; 
    end;
   end




