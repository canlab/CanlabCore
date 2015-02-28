function weig = get_w(sv,method) 

if nargin==1 
    method=1; 
end;
  
weig=[];
for i=1:length(sv.child)  
  weigNew= get_w(sv.child{i});  
  weigNew=weigNew/norm(weigNew);
  weig=[weig;weigNew];
end;


weig=sum(weig);   %% maybe max(weig) could be good too?

