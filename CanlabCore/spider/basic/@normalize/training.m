function [d,a] =  training(a,d)
    
% [results,algorithm] =  training(algorithm,data,loss)
  
  disp(['normalizing with ' get_name(a) '.... '])

  x=get_x(d); l=get_dim(d);

  if a.scale_type==4
     y=get_y(d);
     corr = zeros(size(x,2),1);
     corr = (mean(x(y==1,:))-mean(x(y==-1,:))).^2;
     st   = std(x(y==1,:)).^2+std(x(y==-1,:)).^2;
     a.corr = corr ./ st;
  end
 
 if a.scale_type==1 | a.scale_type==3  | a.scale_type==4
  a.mean_vec  = mean(x);
  a.scale_vec = std(x);
  if a.scale_type==3 %% try to do both scalings
    x = x - ones(l,1)*a.mean_vec; x = x * diag(1./a.scale_vec);
  end
 end
 if a.scale_type==2 | a.scale_type==3
   a.mean_vec2  = mean(x')';
   a.scale_vec2 = std(x')';
 end
     
 
 
 d=test(a,d);
  
  
  
