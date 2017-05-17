function d =  testing(a,d)

 
 r=test(a.child,d); % predict on all bags
 Yest=r{1}.X*0;
 for i=1:a.bags
     Yest=Yest+r{i}.X;
 end
 Yest=Yest/a.bags;  % take mean output over all bags
 
 if a.algorithm.use_signed_output
   Yest=sign(Yest); Yest(Yest==0)=-1;
 end
 
 d=set_x(d,Yest); 
 d=set_name(d,[get_name(d) ' -> ' get_name(a)]); 
  
 
