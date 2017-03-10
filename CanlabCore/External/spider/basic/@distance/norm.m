function D = norm(a,d1,d2,ind1,ind2,kerparam),

  % a bit slow but it works..for now
  
  x=get_x(d2,ind2); y=get_x(d1,ind1);  
  l1=size(x,1); l2=size(y,1); 
  D=zeros(l1,l2);
   
  for i=1:l1
    for j=1:l2
      D(i,j)=norm(x(i,:)-y(j,:),kerparam);
    end
  end  
