
function data_square_plot(a,d)
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
if(nargin <2)
    d=a.Xsv;
end

x=d.X;
y=d.Y;
ax1=min(x); ax2=max(x);
expand = 0.2 * [ax2-ax1];
ax1 = ax1 - expand;
ax2 = ax2 + expand;

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
if 1
    % FeatureLines = [0 -1 1]';  %cheap hack to only get the decision boundary
    colormap('gray')
    % colormap('cool')
    pcolor(gridx, gridy, temp') ;
    [c,h] = contour(gridx, gridy, temp','k') ;

    if(length(h)>1)
        set(h(1),'LineWidth',2);
    end

    i=1;
    while length(h)>i
        set(h(i+1),'LineWidth',2);
        i = i+1;
    end
    % FeatureLines=[-1 1]' ;  % now the SV lines
    [c,h] = contour(gridx, gridy, temp', 'c:') ;
    if ~isempty(h)
        set(h(1),'LineWidth',2,'LineStyle',':');
        i=1;
        while length(h)>i
            set(h(i+1),'LineWidth',2,'LineStyle',':');
            i = i+1;
        end
    end
    shading interp
end


if 1
    pcolor(gridx,gridy,temp');colormap gray;
    [c,h]=contour(gridx,gridy,temp',[-0.1:0.3:1],'k');
    i=1;
    h=get(h,'Children');
    while length(h)>=i
        set(h(i),'LineWidth',3);
        i = i+1;
    end

    if ~isempty(c)
        h=clabel(c,h);
        set(h,'FontSize',16)
    end
else
    mesh(gridx,gridy,temp');colormap gray;
    view(3);
end

shading interp;
colorbar

if(1)
    h=plot(x(c1,1),x(c1,2),'rx'); hold on;
    set(h,'LineWidth',2,'MarkerSize',7);
    h=plot(x(c2,1),x(c2,2),'bo');
    set(h,'LineWidth',2,'MarkerSize',7);
end

sv=find(abs(a.alpha)>1e-7);


%sv=find(abs(a.alpha)>max(abs(a.alpha))/100);
amax=max(abs(a.alpha));

[r]=test(a,a.Xsv);
if(nargin <2)

    if 1
        for i=sv'
            col='co';
            if (r.X(i,:)*r.Y(i,:)<0.9) col='cs';   end;  %% margin error
            h=plot(x(i,1),x(i,2),col);
            alpha=ceil((abs(a.alpha(i))/amax)*4);
            set(h,'LineWidth',alpha,'MarkerSize',8+alpha);
        end
    end
else
    
    Ip=find(y==+1);
    In=find(y==-1);
    
     plot(x(Ip,1),x(Ip,2),'r.','MarkerSize',2)
     plot(x(In,1),x(In,2),'b.','MarkerSize',2)
    
    Z=get_x(a.Xsv);
    plot(Z(:,1),Z(:,2),'go','MarkerSize',2)
end
axis(axis_sz);
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




