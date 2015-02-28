function v =  distortion(a,d)
    
x=d.X; d=test(a,d); Y= d.X;
clusttmp = unique(Y);
K=0;
for i=1:length(clusttmp),
    tmp = find(Y==clusttmp(i));
    K = K + sum((sum((x(tmp,:)-ones(length(tmp),1)*a.mu(clusttmp(i),:)).^2)));     
end;
v=K; 

