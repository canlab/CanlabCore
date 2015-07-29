function data_square_plot(a,granul)
%
% Plots a surface plot of an svm which uses !_2_! dimensional inputs.
% Red decision line is drawn.
%
% Usage :  
%           plot(a)
%
% a -- The algorithm already trained
%      The colormap can be changed afterwards.
%
% Example:
%   d=gen(spiral({'n=1','m=50'}));
%   [r s0]=train(svm(kernel('rbf',1)),d)
%   plot(s0);  


%clf
d=a.child.dat;
x=d.X;
y=d.Y;
ax1=min(x); ax2=max(x);

c1=find(y==1);
c2=find(y==-1);



granul=100;
minX=floor(ax1(1));
maxX=ceil(ax2(1));
minY=floor(ax1(2));
maxY=ceil(ax2(2));
axis_sz=[minX maxX minY maxY];
minX=minX-2; maxX=maxX+2; minY=minY-2; maxY=maxY+2;
gridx=[]; gridy=[];

mx=zeros(granul*granul,2);
for i=1:granul
 for j=1:granul
   mx((i-1)*granul+j,:)= [minX+(i-1)*(maxX-minX)/granul minY+(j-1)*(maxY-minY)/granul ] ;
   gridx(i)=minX+(i-1)*(maxX-minX)/granul;
   gridy(j)=minY+(j-1)*(maxY-minY)/granul;
end
end


temp=zeros(granul,granul);

sflag=a.algorithm.use_signed_output;
a.algorithm.use_signed_output=0;

resX=test(a,data(mx));
resX=sign(get_x(resX)).*sqrt(abs(get_x(resX)));

for i=1:granul
 for j=1:granul
  temp(i,j)= resX( (i-1)*granul+j);
  end
end

surf(gridx, gridy, temp');
view(2);
shading interp





