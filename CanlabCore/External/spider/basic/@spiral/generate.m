function d =  generate(a)


x1=zeros(a.m,1);
y1=zeros(a.m,1);

x2=zeros(a.m,1);
y2=zeros(a.m,1);

step=1/100;
t=step:step:a.n*2*pi;

x1= t.*cos(t) ;
y1= t.*sin(t);
x2= (1+t).*cos(t) ;
y2= (1+t).*sin(t) ;


xfinal1=[];
yfinal1=[];
xfinal2=[];
yfinal2=[];

nn=a.noise;
 for j=1:a.m
     i=ceil(length(t)*rand);
     xfinal1=[xfinal1,x1(i)+nn*rand];
     yfinal1=[yfinal1,y1(i)+nn*rand];
     xfinal2=[xfinal2,x2(i)+nn*rand];
     yfinal2=[yfinal2,y2(i)+nn*rand];
 end

 d=data([get_name(a) ' '] ,[xfinal1' yfinal1';xfinal2' yfinal2'],[ ones(a.m,1);-1*ones(a.m,1)]);