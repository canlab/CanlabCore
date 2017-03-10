function [res,alg] =  training(alg,d,loss_type)
       
  disp(['training ' get_name(alg) '.... '])

  alg.score.child=[]; 
  alg.score.child{1}=alg.child;
  [res, alg.all_algs]=train(alg.score,d);
  %[res,alg.child]=train(cv(alg.child,['folds=' num2str(alg.folds)]),d);
  
  myMean=get_mean; 
  myMean.take_average=1; 
  myMean.loss_type=alg.loss; 
  alg.scores=train(myMean,group(res,'all'));
  alg.scores=group2vec(alg.scores); 
  alg.scores=alg.scores(:,1);
  [alg.best_score,alg.best_index]=min(alg.scores);
  
  global display_tree_array_indexing; % need indexing in tree form 
                                      % to get best of n methods
  isCached=display_tree_array_indexing;  
  display_tree_array_indexing=1;                   
  
  
  % all_algs=alg.all_algs{1};               % remove cv object
  alg.best=alg.all_algs{1}{alg.best_index};         % get best one

%  all_algs=alg.all_algs{1};               % remove cv object
%  alg.best=all_algs{alg.best_index};         % get best one
  
  [res,alg.best]=train(alg.best,d);
  
  
  display_tree_array_indexing=isCached; % revert indexing method back
                                      % to original                   
  



