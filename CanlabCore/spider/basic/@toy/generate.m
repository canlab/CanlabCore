function d =  generate(a)

if a.seed>0 rand('seed',a.seed);   end;
x=(rand(a.l,a.n)-0.5)*2;
w=ones(1,a.n); 
prem=a.n-min(a.relevant_n,a.n);
w(1:prem)=0;
y=sign(x*w');




d=data([get_name(a) ' '] ,x,y);