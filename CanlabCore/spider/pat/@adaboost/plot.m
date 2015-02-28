function plot(a,dd)
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
d=dd;
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


axis(axis_sz);
