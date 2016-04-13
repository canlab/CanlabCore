function plotclustfir(infir,oXc)

colors = {'r' 'g' 'b' 'y' 'c' 'm' 'k'};
if max(oXc)>length(colors)
    error('not enough colors!');
end

for c=1:size(infir,1);
    for s=1:size(infir,2);
        fir(c,:,:,s)=infir{c,s};
    end
end

mfir=squeeze(mean(fir,4));
figure;
for c=1:max(oXc);
leg{c}=['clust',num2str(c)];
end

for n=1:size(fir,1);        % for each condition
    subplot(1,size(mfir,1),n);
    for c=1:max(oXc);
    hold on;
    plotfir=squeeze(mean(mfir(n,:,oXc==c),3));
    plotfir=squeeze(mfir(n,:,oXc==c));    
    plot(plotfir,'color',colors{c},'LineWidth',3)
    end
    legend(leg)
end


