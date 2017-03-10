function [d,a] =  training(a,d)
  
  if a.algorithm.verbosity>0
    disp(['training ' get_name(a) '.... '])
  end
  
  a.mean=mean(d.X);
      
  if strcmp(a.assume,'equal_cov')
    a.cov=mean(var(d.X));
  else
    if strcmp(a.assume,'diag_cov')
      a.cov=var(d.X);
    else
      a.cov=cov(d.X);  
    end
  end
  
  if ~a.algorithm.do_not_evaluate_training_error
    d=test(a,d);
  end
  
  
