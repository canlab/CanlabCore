
% -----------------------------------------------------------------------
% plot output
% -----------------------------------------------------------------------
%%%%%%%%%%bar graph
% scnsize = get(0,'ScreenSize');
% myposition = [50 50 scnsize(3)-100 scnsize(4)/2];
% myposition(3) = min(myposition(3), 1200);
% myposition(4) = min(myposition(4), 500);
% f1 = figure('position', myposition,'color','white'); set(gca,'FontSize',14),hold on
% subplot(1,3,1); set(gca,'FontSize',12);


% Programmers notes:
% 8/21/2017 - Stephan
% changed xticklabels to auto-labeling. The previous version had ticks for
% every cluster, which got unreadable for k>15. 


f1 = create_figure('TestClust_Quality',1,3);


%%%------------------------- panel 1 ---------------------------------%%%
cmcq = cmcq(clust); %return to original format
cstdval = cstdval(clust);
% bestc = bestc - clust(1)+1;
xpos  = 2:length(clust)+1;
hold on;
plot(xpos,cmcq,'ko-','LineWidth',3,'MarkerFaceColor','k');
if exist('bestmcq', 'var')
    plot(bestc,bestmcq,'ks','LineWidth',3,'MarkerFaceColor',[.5 .5 .5],'MarkerSize',13);
end

plot(xpos,permuted_quality,'o-','Color',[.5 .5 .5],'MarkerFaceColor','k','LineWidth',2);
plot(xpos,perm95,'--','Color',[.5 .5 .5],'LineWidth',2);
plot(xpos,perm05,'--','Color',[.5 .5 .5],'LineWidth',2);

legend({'Real-data solutions' 'Chosen solution' 'Permuted solutions' '95% confidence'},'location','best');
legend boxoff;
set(gca,'XLim',[.5 length(clust)+1.5])

%if ~isnan(coq)
%    set(gca,'Ylim',[0 max(coq)+max(cstd_poq)+.5]);
%end

title('Cluster quality')
xlabel('Number of clusters in solution')
ylabel('Mean silhouette value')
%%%------------------------- panel 1 ---------------------------------%%%


%%%------------------------- panel 2 ---------------------------------%%%
subplot(1,3,2); %set(gca,'FontSize',12);
hold on;
plot(xpos,cstdval,'ko-','LineWidth',3,'MarkerFaceColor','k')
plot(bestc,beststdval,'ks','LineWidth',3,'MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)

%plot([bestc bestc],[min(cstdval) max(cstdval)-.01],'k:','Linewidth',2);

legend({'Real-data solutions' 'Chosen solution'},'location','best');
legend boxoff;
set(gca,'XLim',[.5 length(clust)+1.5])

xlabel('Number of clusters in solution')
ylabel('Improvement over permuted data (s.d.)')
%%%------------------------- panel 2 ---------------------------------%%%


%%%------------------------- panel 3 ---------------------------------%%%
%%%%%%%%%%%%plot best histogram
subplot(1,3,3); %set(gca,'FontSize',12); hold on;
%for c=1:length(clust);

[h,x] = hist(bestpmcq,min(10,round(nperm ./ 40)));      %plot hist of best distribution
hh = bar(x,h);
set(hh,'FaceColor',[.3 .3 .3],'EdgeColor',[.3 .3 .3])

if exist('bestmcq', 'var')
    x = [bestmcq bestmcq];      %0:max(h);
    y = [0 max(h)];                  %ones(1,max(h)+1); %   *coq(best);

    plot(x,y,'color','k','LineWidth',3)
end

%set(gca,'XLim',min(h) - .1
title(['Permuted distribution for ',num2str(bestc),' clusters'])
ylabel('Mean silhouette value')
ylabel('Frequency')
%%%------------------------- panel 3 ---------------------------------%%%