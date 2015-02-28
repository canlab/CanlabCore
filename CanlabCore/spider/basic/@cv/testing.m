function [retRes] = testing(algo,dat)

% Does cross validation and returns the results in a data object
% plus the updated group (after training models)

     
  if algo.store_all==0
   error('cannot test cv objects as the models have not been stored, use store_all=1 hyperparameter on cv instead');
  end
  retRes=[];
  for i=1:length(algo.child)
    res=test(algo.child{i},dat);
    res=make_cell(res);
    if ~isempty(retRes) 
        retRes={retRes{1:length(retRes)} res{1:length(res)}}; 
    else 
        retRes=res; 
    end
  end  
 
  if length(retRes)==1 
      retRes=retRes{1}; 
  else 
      retRes=group(retRes); 
  end;
  


