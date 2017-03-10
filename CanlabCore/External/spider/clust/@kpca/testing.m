function d =  testing( a, d)

  K = calc( a.child, d, a.dat)';   %% between old examples and test examples
  if 0
    [Kt] = calc(a.child,a.dat,a.dat); %%restore old uncentered kernel
  else
    Kt=a.Kt;
  end
  
  if a.center_data   %center the test examples in feature space
    [n,m] = size(K);
    Om = ones(m)/m;
    On = ones(n,m)/m;
    I = eye( m); 
    %brackets are for better numerical condition
    K = ((K - On*Kt)*(I-Om))';

    test_features = K'*a.e_vec;
else
    test_features = K*a.e_vec;
end
  d=set_x(d,test_features); 
  d.name=[get_name(d) ' -> ' get_name(a)];
  
