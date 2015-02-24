function x=makebinary(y);
y=y(:)';
x=zeros(max(y),length(y));

for n=1:max(y);
    x(n,:)=y==n;
end
    
    