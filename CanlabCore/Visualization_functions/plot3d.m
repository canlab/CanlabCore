function plot3d(X,names,clus,linemat,cordir)

%%% define scaling %%%%%%%
numel=size(X,1);
a=X(:,1);b=X(:,2);c=X(:,3);
axa=ones(numel,1)*mean(a);
axb=ones(numel,1)*mean(b);
axc=ones(numel,1)*mean(c);
axa1=min(a):(max(a)-min(a))/(numel-1):max(a);
axb1=min(b):(max(b)-min(b))/(numel-1):max(b);
axc1=min(c):(max(c)-min(c))/(numel-1):max(c);

%%% plotting defaults %%%%
f1=figure;
hold on;
set(gca,'FontSize',10);
colors = {'ko' 'ro' 'go' 'bo' 'yo' 'co' 'mo' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
while max(clus) > length(colors), colors = [colors colors];,end

%draw central axes & scale
plot3(axa,axb,axc1,'color','k','linewidth',3);
plot3(axa,axb1,axc,'color','k','linewidth',3);
plot3(axa1,axb,axc,'color','k','linewidth',3);
axis([min(a) max(a) min(b) max(b) min(c) max(c)]);

%%% plot lines %%%%%%
for r1=1:size(linemat);
    for r2=r1+1:size(linemat);
        if cordir(r1,r2)==1;
            col=[0.7 0.7 0.7];
        else
            col=[0 0 0.7];
        end
        if linemat(r1,r2)==1;
           plot3([X(r1,1) X(r2,1)],[X(r1,2) X(r2,2)],[X(r1,3) X(r2,3)],'color',col,'linewidth',2);
        end
    end
end

%%%% plot points %%%%
for j = 1:size(X,1);
    plot3(a(j),b(j),c(j),colors{clus(j)},'MarkerFaceColor',colors{clus(j)}(1),'Markersize',8);
end

if ~isempty(names);
for i = 1:size(X,1)
    text(X(i,1)+.025*[max(X(:,1))-min(X(:,1))],X(i,2)+.0,X(i,3)+.0,names{i},'Color','k','FontSize',12,'FontWeight','bold');
end
end

rotate=0;
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

figure;hold off;
subplot(2,2,1);
mdsfig(X(:,[1 2]),names,clus,linemat);
subplot(2,2,2);
mdsfig(X(:,[1 3]),names,clus,linemat);    
subplot(2,2,3);
mdsfig(X(:,[2 3]),names,clus,linemat);    



