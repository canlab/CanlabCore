function data_square_plot(a,d)
% clf;
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





% clf


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

hold on;

%surf(1:granul,1:granul,temp'-1000);







%clf

FeatureLines = [0 -1000 1000]';  %cheap hack to only get the decision boundary
h = surf(gridx,gridy,temp') ;


shading interp







%surf(gridx,gridy,temp');

%view(2);

%shading interp; 

%[c,h]=contour(gridx,gridy,temp',[0 0],'k');

%   clabel(c,h);

%colorbar



%h=plot(x(c1,1),x(c1,2),x(c1,2)*0+100,'rx'); hold on;

%set(h,'LineWidth',2,'MarkerSize',7);

%h=plot(x(c2,1),x(c2,2),x(c2,2)*0+1000,'bo');

%set(h,'LineWidth',2,'MarkerSize',7);


 



axis off;

% 

% if( exist('dat'))

%     for i=1:length(dat.X(:,1))

%         if dat.Y(i,1) > 0

%             plot3( granul+granul*(minX+dat.X(i,1))/(max_x-minX) ,granul+granul*(minY+dat.X(i,2))/(maxY-minY),-100,'resX');

%         else

%             plot3( granul+granul*(minX+dat.X(i,1))/(max_x-minX) ,granul+granul*(minY+dat.X(i,2))/(maxY-minY),-100,'go');

%         end

%     end

% end









