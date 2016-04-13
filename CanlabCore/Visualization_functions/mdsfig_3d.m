function mdsfig_3d(X,names,clus,linemat,cordir)
% ::
%
%    mdsfig_3d(X,names,clus,linemat,cordir)
%
% 3-D MDS figure
%
% Subfunction of mdsfig used if 3+ dims are available

numel=size(X,1); %%% define scaling %%%%%%%
a=X(:,1);b=X(:,2);c=X(:,3);
axa=ones(numel,1)*mean(a);
axb=ones(numel,1)*mean(b);
axc=ones(numel,1)*mean(c);
axa1=min(a):(max(a)-min(a))/(numel-1):max(a);
axb1=min(b):(max(b)-min(b))/(numel-1):max(b);
axc1=min(c):(max(c)-min(c))/(numel-1):max(c);

%%% plotting defaults %%%%
f1=figure('Color','w');
hold on;
set(gca,'FontSize',12);
colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
while max(clus) > length(colors), colors = [colors colors];,end

%draw central axes & scale
plot3(axa,axb,axc1,'color','k','linewidth',2);
plot3(axa,axb1,axc,'color','k','linewidth',2);
plot3(axa1,axb,axc,'color','k','linewidth',2);
axis([min(a) max(a) min(b) max(b) min(c) max(c)]);

%%% plot lines %%%%%%

drawlines(X,linemat,0,{'Positive' 'Negative'});


%%%% plot points %%%%
for j = 1:size(X,1);
    plot3(a(j),b(j),c(j),colors{clus(j)},'MarkerFaceColor',colors{clus(j)}(1),'Markersize',8);
end

if ~isempty(names);
for i = 1:size(X,1)
    text(X(i,1)+.025*[max(X(:,1))-min(X(:,1))],X(i,2)+.0,X(i,3)+.0,names{i},'Color','k','FontSize',12,'FontWeight','bold');
end
end

rotate=1;
if rotate
% rotate until closed
for v=1:45;
    if ishandle(f1);
        view(v,35);
        drawnow;
        pause (0.2)
        if v==44
            v==1;
            continue
        end
    else
        break
    end
end
end


% Second Figure -- different dims

tor_fig(2,2);hold off;
subplot(2,2,1);
mdsfig(X(:,[1 2]),names,clus,linemat);
subplot(2,2,2);
mdsfig(X(:,[1 3]),names,clus,linemat);    
subplot(2,2,3);
mdsfig(X(:,[2 3]),names,clus,linemat);    

subplot(2,2,4);
if size(X,2) > 3, 
    mdsfig(X(:,[1 4]),names,clus,linemat);  
else
    axis off
end


return







function [hhp,hhn] = drawlines(pc,sigmat,sigcol,legmat,varargin)
    
% sigcol is a scalar between zero and one, lower is 'more salient'
hold on

% linew
lw = 3;

color(1,:) = [0 0 0];   % first line color, 'positive'
color(2,:) = [0 .7 1];   % second line color, 'negative'
style{1} = '-';         % first line style
style{2} = '-';         % second line style

if length(varargin) > 0, color = varargin{1}; end
if length(varargin) > 1, style = varargin{2}; sigcol = 0; end

hhp=[];hhn=[];
for i=1:size(pc,1);
     for j=1:size(pc,1);
         if sigmat(i,j) > 0
             hhp = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
             set(hhp,'Color',color(1,:) + [1 1 1] * abs(sigcol),'LineStyle',style{1},'LineWidth',lw - (lw*sigcol)-1)
         elseif sigmat(i,j) < 0
             hhn = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
             set(hhn,'Color',color(2,:) + abs([sigcol sigcol 0]),'LineStyle',style{2},'LineWidth',lw - (lw*sigcol)-1)
         end
     end
end
legend ([hhp hhn],legmat);
return


