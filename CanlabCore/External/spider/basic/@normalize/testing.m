function d =  testing(a,d)

 x=get_x(d); [l n]=get_dim(d);
 
 if a.scale_type==1 | a.scale_type==3 | a.scale_type==4
   x = x - ones(l,1)*a.mean_vec;
   x = x * diag(1./a.scale_vec);
 end
 if a.scale_type==2 | a.scale_type==3 
   x = x - a.mean_vec2*ones(1,n);
   x = (x' * diag(1./a.scale_vec2))';
 end
 
 if a.sigmoid==1
   x=tanh(x*a.sig);
 end
 
 if a.scale_type==4
     x = x * diag(a.corr);
 end     
   
 d=set_x(d,x);
 d=set_name(d,[get_name(d) ' -> ' get_name(a)]);
 
 

