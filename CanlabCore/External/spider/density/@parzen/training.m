function [d,a] =  training(a,d)
  
  if a.algorithm.verbosity>0
    disp(['training ' get_name(a) '.... '])
  end

  a.child.do_not_evaluate_training_error=1;
  [r a.child]=train(a.child,d); %% make kernel store data
  
  if ~a.algorithm.do_not_evaluate_training_error
    d=test(a,d);
  end
  
  
