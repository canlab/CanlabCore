function d =  testing(a,d)
  
 d=test(a.child,d);
 d.X=mean(d.X)';
 
 
 d=set_name(d,[get_name(d) ' -> ' get_name(a)]); 
  
 
