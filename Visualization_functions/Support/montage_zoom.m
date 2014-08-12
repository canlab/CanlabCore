function montage_zoom
% montage_zoom, no arguments
% in a montage created with montage_clusters,
% zooms in on current axis
%
% tor wager

todel = get(gcf,'Children');
todel = todel(~(todel == gca));

%tmp = find(todel == gca);
%try,delete(todel(1:tmp-2)),catch,end

set(todel,'Position',[0.01 0.01 0.02 0.02])

%todel = get(gcf,'Children');
%todel = todel(~(todel == gca));
%delete(todel)

set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])

h = findobj(gca,'Type','text');
set(h,'FontSize',24)