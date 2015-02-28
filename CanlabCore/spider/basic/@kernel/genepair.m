function K = poly(kern,dat1,dat2,ind1,ind2,kerParam),
 
% K = poly(d1,d2,ind1,ind2,param), compute the kernel matrix between d1 and d2
%  for a polynomial kernel (<x,z>+1)^param where x is from d1 and z from d2
  
K=genepair_calc(get_x(dat2,ind2),get_x(dat1,ind1));

 
