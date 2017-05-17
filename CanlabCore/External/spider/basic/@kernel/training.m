
function [d,a] = training(a,d)

  a.dat=d;

  if ~a.algorithm.do_not_evaluate_training_error
    d=set_x(d,feval('calc',a,d)); 
  end
