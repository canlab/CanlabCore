function dat =  testing(a,dat)
 
 feat=a.feat; if isempty(feat) feat=length(a.rank);   end;
 
 
 gi = a.feature_group;
 if isempty(gi)
	[ans, nfeatures] = get_dim(dat);
	gi = 1:nfeatures;
 end
 select = find(ismember(gi, a.rank(1:feat)));
 
 
 dat=get(dat,[],select);  % perform the feature selection
 dat=set_name(dat,[get_name(dat) ' -> ' get_name(a)]);
 
 if a.output_rank==0  %% output selected features, not label estimates
   dat=test(a.child,dat);  % train underlying algorithm 
 end
 