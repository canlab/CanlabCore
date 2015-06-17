function h = plotDelta(delta)
% h = plotDelta(delta)
% 
% A plot of a random ER design:
% 
% stim = [ones(20,1); 2*ones(20,1); 3*ones(20,1); 4*ones(20,1); 0*ones(100,1)]; stim =getRandom(stim);
% [X,d] = getPredictors(stim,spm_hrf(1)./max(spm_hrf(1)));
% tor_fig; plotDelta(d);
% axis on; set(gca,'YColor',[1 1 1],'FontSize',16);
% xlabel('Time (s)')
%
% tor_fig; hrf = spm_hrf(.1)./max(spm_hrf(.1)); 
% x = 0:.1:32; plot(x(1:length(hrf)),hrf,'k','LineWidth',2); 
% axis on; set(gca,'YColor',[1 1 1],'FontSize',16);
% xlabel('Time (s)')
%
% tor_fig;
% for i = 1:L
% subplot(L,1,i)
% hold on
% plot(X(:,i),'k','LineWidth',2);
% axis off,set(gca,'YLim',[-.5 3])
% end
% axis on; set(gca,'YColor',[1 1 1],'FontSize',16);
% xlabel('Time (s)')


if ~iscell(delta)
    for i = 1:size(delta,2), d{i} = delta(:,i);, end
    delta = d;
end

L = length(delta);

for i = 1:L                                
    subplot(L,1,i)                               
    hold on

    d = find(delta{i});
    z = zeros(size(d));

    plot([d d]',[z z+1]','k','LineWidth',2);

    set(gca,'YLim',[-.5 1])
    axis off
end

return
