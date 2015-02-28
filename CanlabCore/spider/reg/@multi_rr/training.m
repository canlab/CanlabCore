function [results,a] =  training(a,d)  
         
  % [results,algorithm] =  training(algorithm,data,loss)  
  
  disp(['training ' get_name(a) '.... '])  
    
  if ~isempty(a.indices) & ~a.use_kernels  
    d=get(d,a.indices);  
  end  
 
  %d2=get(d,a.indices);K =get_kernel(a.child,d,d2);  %% OLD CODE
  
  a.Xsv=d;
  if a.use_kernels
   [K,a.child] =train(a.child,d); K=K.X;  
  else
   K=d.X;
  end   
  
  if(a.use_b & ~a.use_kernels) K=[K; ones(1,size(K,2))]; end;  

  if a.use_kernels 
    Kinv=K; Kinv=Kinv+eye(length(Kinv))*a.ridge;  
    
    if (rcond(Kinv) > 1e-15)
        Kinv=inv(Kinv); 
    else
        disp(['Bad conditioned matrix cond(K + r*I) = ' num2str(rcond(Kinv)) '! Taking pseudo inverse...']);
        Kinv=pinv(Kinv);
    end
    
    Y = get_y(d);
    a.alpha=Kinv*Y;  %% calculate alphas   
  else 
    Kinv=K'*K; Kinv=Kinv+eye(length(Kinv))*a.ridge;  
    
    if (rcond(Kinv) > 1e-15)
        Kinv=inv(Kinv); 
    else
        disp(['Bad conditioned matrix cond(K + r*I) = ' num2str(rcond(Kinv)) '! Taking pseudo inverse...']);
        Kinv=pinv(Kinv);
    end
    
    Y = get_y(d);
    a.alpha=(Kinv*K'*Y);  %% calculate alphas   
  end
 
    
  results=test(a,d);  
    
  
  
  
  
  
  
