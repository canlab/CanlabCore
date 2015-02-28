function [r,a] =  training(a,d)
  
  if a.algorithm.verbosity>0
    disp(['training ' get_name(a) '.... '])
  end
    
  c=size(get_y(d),2); % #classes
  n=size(get_x(d),1);
  W=ones(n,1)/(c*n);
  for k=1:a.kmax
  	cW = cumsum(W);
	randnum = rand(n,1);
	m=n;
	indices = zeros(m,1);	
   	for i = 1:m
      		%Find which bin the random number falls into
      		loc = max(find(randnum(i) > cW));		
      		if isempty(loc)
         		indices(i) = 1;
      		else
         		indices(i) = loc;
      		end
   	end	
  	[r,a.child{k}] =  train(a.child{k},get(d,indices));
	E = get_y(loss(r,'weighted_class_loss',W(indices,:)));  
	if(E > 0.5) 		
		k=k-1;	
		break;
	end   	
	if(E == 0)
		E=0.001;
	end
	a.alpha(k) = log((1-E)/E);	
	W  = W.*exp(-a.alpha(k)*max((r.X~=r.Y)')');
   	W  = W./sum(W);	
  end
  a.alpha=a.alpha./sum(abs(a.alpha(1:k)));  
  a.kmax=k;
  
  if ~a.algorithm.do_not_evaluate_training_error
    r=test(a,d);
  else
    r=d;
    r=set_x(r,get_y(r));
  end
  
  
