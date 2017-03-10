function K = rbf_of_dist(kern,dat1,dat2,ind1,ind2,kerParam),
  %
  % introduce nonlinearity by performing POLYNOMIAL on input kernel matrix
  %
  K=get_x(dat1,ind1);
  K = exp(-K/(2*kerParam^2));
  
