function d =  testing(a,d)
  

dist=calc(a.child,data(a.mu),d);  
  
if a.algorithm.use_signed_output
  [dummy,Yest] = min(dist');
end
  
d=set_x(d,Yest');  d.name=[get_name(d) ' -> ' get_name(a)]; 
