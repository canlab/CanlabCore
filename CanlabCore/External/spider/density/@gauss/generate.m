
function d =  generate(a)

x=randn(a.l,length(a.mean));
n=length(a.mean);

if size(a.cov,2)==1
  cov=eye(n)*a.cov;   %% make cov matrix with only diagonal elems (all same)
end;
if size(a.cov,1)==1
  cov=diag(a.cov);   %% make cov matrix with only diagonal elems (different) 
end
  

for i=1:a.l
  if length(a.cov)==1
    x(i,:)= (x(i,:)*sqrt(a.cov)) + a.mean; continue;
  else
    if size(a.cov,1)==1
      x(i,:)= (x(i,:).* sqrt(a.cov)) + a.mean; continue;
    else
      x(i,:)= x(i,:) * a.cov^0.5 + a.mean;
    end
  end 
end

d=data(['gauss' ' l=' num2str(a.l)],x,[]); 
