function [r,a] =  training(a,d) 
  
% [results,algorithm] =  train(algorithm,data)

  store=[];
  a.store_all=1; %%%% train only useful for testing if you store models
		 
  if isa(a.child,'cell') a.child=a.child{1}; end; 
  a.child.algorithm.trained=0;
   
  if ~isa(a.values,'cell') val{1}=a.values; else val=a.values; end
  if ~isa(a.param,'cell') p{1}=a.param; else p=a.param; end
   
  for i=1:length(val) sz(i)=length(val{i}); end;
  tot=prod(sz); %% find permutations of hyperparameters 
  r=[];  
  for i=1:tot %% all permutations
    vars=num2choice(a,i,sz); %% hyperparams for this iteration      
    for j=1:length(p)      %% set hypers       
      v=(val{j}(vars(j)));  
      
      f=a.child; ll=length(f); %length(get(f,'child'));
      if ll==1
	if strcmp(a.algorithm.name,'param') %% need this explicit
                %for param object as matlab doesnt use subasgn
                %because we're inside the param object here in the code
	  s=[]; s.type='.'; s.subs=p{j};
	  a.child=subsasgn(f,s,v);
	else
	  eval([ 'a.child.' p{j} '=v;']);
	end
      else 
	for z=1:ll
	  eval([ 'a.child{z}.' p{j} '=v;']); 
	end
      end
    
      %% set whether to force/ not force training for that object
      e=p{j}; t=max(find(e=='.')); f=a.force_no_train;
      eval(['a.child.' e(1:t) 'algorithm.do_not_retrain=f;']);
      %% if not the first run then re-use/not re-use what has been done
      %% before
     if ~isempty(r),
        f=a.force_use_prev_train;
        eval(['a.child.' e(1:t) 'algorithm.do_not_retrain=f;']);
     end; 
    end 
    
    if isa(d,'group') & strcmp(d.group,'one_for_each') & length(d.child)==length(a.child)
      [res, a.child]=train(a.child,d.child{i}); %% special 'one_for_each' group
    else
      [res, a.child]=train(a.child,d);
    end
    
    
    if a.store_all  %% store results
      store{i}=a.child;
    end;
    r{i}=res;
  end
  
  if a.store_all  %% store results
   a.child=store; 
  end

  
  %% reset the values of use_prev_train  
  if a.store_all,
    for i=1:length(a.child),
       for j=1:length(p),
           e=p{j}; t=max(find(e=='.'));
	   eval(['a.child{i}.' e(1:t) 'algorithm.use_prev_train=0;']);	   
       end;
    end;
  else
      for j=1:length(p),
           e=p{j}; t=max(find(e=='.'));
	   eval(['a.child.' e(1:t) 'algorithm.use_prev_train=0;']);	   
      end;
  end;

  if length(r)==1 r=r{1}; else r=group(r); end;
  if isa(r,'group') r.group=a.group; end;
  


  %a=group(a.child);  %% could just do this, actually...




