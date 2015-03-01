
function [ret,algo] = training(algo,dat)
  
  ret=feval(algo.type,algo,dat); 
