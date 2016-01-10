function hewma_cluster_plot(cl,i)

ti = cl(i).timeseries;

figure; fa = gca; set(fa,'FontSize',16)
plot(ti,'k','LineWidth',2);
hold on; fill_around_line(ti,cl(i).timeseries_ste,'k');

% title and cp line
p = 'hewma_t.img';
tmp = tor_extract_rois(p,cl(i)); t = tmp.timeseries;
p = 'hewma_z.img';
tmp = tor_extract_rois(p,cl(i)); z = tmp.timeseries;

p = 'hewma_cp.img';
tmp = tor_extract_rois(p,cl(i)); cp = tmp.timeseries;
p = 'hewma_p.img';
tmp = tor_extract_rois(p,cl(i)); pval = tmp.timeseries;

str = sprintf('Cl %3.0f, t=%3.2f, zcorr=%3.2f,p=%3.4f',i,t,z,pval);

axes(fa)
title(str)
plot([cp cp],get(gca,'YLim'),'k--');


% slice image

if mean(ti(1:10)) > mean(ti(end-9:end))
    a = axes('Position',[.65 .65 .25 .25]);
else
    a = axes('Position',[.1 .65 .25 .25]);
end
drawnow

montage_clusters_maxslice([],cl(i),{'r'});
drawnow


return
