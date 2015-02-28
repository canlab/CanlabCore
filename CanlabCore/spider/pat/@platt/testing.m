function d =  testing(a,d)

r=test(a.child,d);
out = get_x(r);
out = 1./(1+exp(a.A*out + a.B));

if a.algorithm.use_signed_output==1
   Yest=sign(out-0.5);
else
    Yest = out;
end

 
d=set_x(d,Yest); 
d=set_name(d,[get_name(d) ' -> ' get_name(a)]); 
  
