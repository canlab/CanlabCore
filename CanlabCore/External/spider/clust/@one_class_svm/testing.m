function d =  testing(a,d)
  
 Kt=get_kernel(a.child,d,a.Xsv);
 Yt=get_y(d); 
 
 Yest=((a.alpha'* Kt)+a.b0)';
  
 if a.algorithm.use_signed_output==1
   Yest=sign(Yest);
 end 

 d =set_x(d,Yest); d.name=[get_name(d) ' -> ' get_name(a)]; 
 
 
 
