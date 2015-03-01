function K = poly_of_ker(kern,dat1,dat2,ind1,ind2,kerParam),

  % introduce nonlinearity by performing POLYNOMIAL on input kernel matrix
  
  K=get_x(dat1,ind1);
  K= (K+1).^kerParam;
  
  
