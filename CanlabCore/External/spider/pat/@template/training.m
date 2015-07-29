function [d,a] =  training(a,d)
  
  if a.algorithm.verbosity>0
    disp(['training ' get_name(a) '.... '])
  end
  
  X=get_x(d); K=X*X'; Y=get_y(d); K=K.*(Y*Y');
  [alpha,b] = quadsolve(K,-ones(size(K,1),1),Y',0,a.C); 
  alpha= alpha .* Y;
  
  a.alpha=alpha;
  a.b0=-b;
  a.Xsv=d;

  if ~a.algorithm.do_not_evaluate_training_error
    d=test(a,d);
  end
  
  
