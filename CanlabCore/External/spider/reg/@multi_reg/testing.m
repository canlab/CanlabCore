function d =  testing(a,d)
 

for i=1:a.Q,
  Yest(:,i) = get_x(test(a.child{i},d));    
end;
  
d.name=[get_name(d) ' -> ' get_name(a)]; d.X=Yest;
