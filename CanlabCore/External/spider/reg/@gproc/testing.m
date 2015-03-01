function dat =  testing(algo,dat)
 pred = get_x(dat); 
 
%  if algo.pred_var
%  	[yEst s2] = gpS00(algo.H, {algo.input, algo.target, pred});
% 	dat=set_x(dat,s2); 
%  else
  	[yEst s2] = gpS00(algo.H, {algo.input, algo.target, pred});
	dat=set_x(dat,[yEst,s2]); 
%  end 
 
 dat=set_name(dat,[get_name(dat) ' -> ' get_name(algo)]); 
  
