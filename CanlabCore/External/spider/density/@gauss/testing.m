function d =  testing(a,d)
   
 if isempty(d) %% generate data according to model instead
   d=generate(a); return;
 end
  
 [l n k]=get_dim(d);

 %% --- now convert to a cov matrix if only storing diagonal elements
 
 cov=a.cov; 
 if size(a.cov,2)==1
   cov=eye(n)*a.cov;   %% make cov matrix with only diagonal elems (all same)
 end;
 if size(a.cov,1)==1
   cov=diag(a.cov);   %% make cov matrix with only diagonal elems (different) 
 end
 invcov=inv(cov + eye(n)*1e-10 );
 
 x=d.X; Yest=zeros(size(x,1),1);
 n=1/( (2*pi)^(n/2) * norm(a.cov)^0.5 );
  
 for i=1:size(x,1)
   xs= (x(i,:) - a.mean);
   xs= xs*invcov*xs';  
   Yest(i)=n*det(cov)^0.5 * exp(-0.5*xs);
 end

 d=set_x(d,Yest); 
 d=set_name(d,[get_name(d) ' -> ' get_name(a)]); 
  
 
