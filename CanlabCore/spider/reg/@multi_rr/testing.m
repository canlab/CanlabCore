function d =  testing(a,d)  
    
  %K = get_kernel(a.child,d,a.X); %% OLD CODE 

  if a.use_kernels
    K = test(a.child, d);  K=K.X; 
  else
    K=d.X;
  end 

  if(a.use_b & ~a.use_kernels) K=[K ; ones(1,size(K,2))]; end;       
  Yest = K'*a.alpha;     
 
  if a.algorithm.use_signed_output==1  
    Yest=sign(Yest);  
  end  
  d=set_x(d,Yest);    
  d=set_name(d,[get_name(d) ' -> ' get_name(a)]);   
   
   
  
  
   
  
  
  
  
