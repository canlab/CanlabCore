function K = rbf(kern,dat1,dat2,ind1,ind2,kerParam),

% K = rbf(d1,d2,ind1,ind2,param), compute the kernel 
%     matrix between d1 and d2
% for a rbf kernel exp(-||x-z||^2/(2*param^2)) 
%     where x is from d1 and z from d2

  K=get_x(dat2,ind2)*get_x(dat1,ind1)';  
  kernTemp=kernel;
  Kdn = get_norm(kernTemp,dat1,ind1).^2; 
  Kn = get_norm(kernTemp,dat2,ind2).^2;  
  K = ones(length(Kn),1)*Kdn' + Kn*ones(1,length(Kdn)) - 2*K;
  K = exp(-K/(2*kerParam^2));
