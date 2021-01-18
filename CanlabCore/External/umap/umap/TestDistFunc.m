function D2=TestDistFunc(Z1,ZJ)
[R,C]=size(ZJ);
D2=zeros(R,1);
for r=1:R
    for c=1:C
        D2(r)=Z1(c)-ZJ(r,c);
    end
end
end
