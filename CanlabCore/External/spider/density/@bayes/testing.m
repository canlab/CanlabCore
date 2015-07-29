function d =  testing(a,d)
  
 if isempty(d) %% generate data according to model instead
   d=generate(a); return;
 end
 
 [l n k k2]=get_dim(d);
 Yest=zeros(l,k2);
 
 for i=1:length(a.child) 
   Yest(:,i)=a.prior(i)*get_x(test(a.child{i},d));
 end

 if a.algorithm.use_signed_output
   [a1 a2]=max(Yest');
   Yest=convert_mc(a2',k2);
   if size(Yest,2)==2 Yest=Yest(:,1); end;
 end
 
 d=set_x(d,Yest); 
 d=set_name(d,[get_name(d) ' -> ' get_name(a)]); 
  
 




