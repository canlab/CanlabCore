function D = linear(k,d1,d2,ind1,ind2,kerparam),

  D=get_x(d2,ind2)*get_x(d1,ind1)';  
  d=distance;
  Ddn = get_norm(d,d1,ind1).^2; 
  Dn = get_norm(d,d2,ind2).^2;  
  D = ones(length(Dn),1)*Ddn' + Dn*ones(1,length(Ddn)) - 2*D;
  D = sqrt(D);