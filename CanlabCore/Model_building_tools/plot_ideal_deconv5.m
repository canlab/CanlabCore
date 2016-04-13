function plot_ideal_deconv5(rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,INFO,hrf,snr,TR,varargin)
% :Usage:
% ::
%
%     function plot_ideal_deconv5(rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,INFO,hrf,snr,TR,[truer])
%
% delta should be matrix of column vectors
%
% optional: truer, vector of "true" responses for each trial type
%
% ..
%    Tor Wager
% ..

msdblim = max(max(abs(msdb)));
bmlim = max(max(abs(biasmean))); bmlim = [-bmlim bmlim];

if isfield(INFO,'ideal_data') & isfield(INFO,'delta')
    figure('Color','w'), myleg{1} = 'Ideal response';
    h(1) = plot(INFO.ideal_data,'LineWidth',2);, hold on;
    mycols = {'r' 'g' 'y' 'c' 'm'}; mycols = repmat(mycols,1,10);
    delta = INFO.delta;
    for i = 1:size(delta,2)
        tmp1 = find(delta(:,i));
        hh = plot([tmp1 tmp1]',[-1*ones(size(tmp1)) zeros(size(tmp1))]',mycols{i},'LineWidth',2);
        h(i+1) = hh(1);
        myleg{i+1} = ['Trial type ' num2str(i) ' onsets'];
    end
    set(gca,'FontSize',18)
    legend(myleg)
    set(gcf,'Position',[19         606        1561         360])
    
    %tmp = conv(hrf,delta(:,1));figure;plot(tmp),title('ttype1')
    %tmp = conv(hrf,sum(delta,2));figure;plot(tmp),title('alltypes')
end
%figure;plot([find(tmp);find(tmp)],[zeros(size(find(tmp)),1); ones(size(find(tmp)),1)])

if length(varargin) > 0
    truer = varargin{1};
else
    truer = ones(1,size(INFO.delta,2));
end


figure; subplot 221; hold on; set(gcf,'Color','w');
if ~isempty(INFO.autocorrelation) 
    figure('Color','w'); set(gca,'FontSize',18)
    plot((1:length(INFO.autocorrelation))./(length(INFO.autocorrelation) / length(hrf)),INFO.autocorrelation,'k-','LineWidth',2)
end
title('Noise autocorrelation function')

% old plot of RMSD
% -----------------
%plot(snr,rmsd,'ko-','MarkerFaceColor','k','MarkerSize',3,'LineWidth',2)
%plot(snr,msdstd(:,1),'k--')
%plot(snr,msdstd(:,2),'k--')
%legend({'Mean population value' '95% CI for one regression'})
%ylabel('Root Mean Squared Deviation')
%xlabel('Signal to Noise Ratio')

%title(['Est. - True response: ' INFO.condition_tested],'FontSize',14)
 
subplot 222; hold on;
x = (1:length(hrf)) * TR;
plot(x, hrf.*truer(INFO.ttype),'k','LineWidth',2)
myleg = {'HRF'};
legend(myleg)
title(['Est. vs. True response: ' INFO.condition_tested],'FontSize',14)
xlabel('Time from trial onset')

mycolors = {'r' 'b' 'g' 'c' 'm'};
for i = 1:min(5,length(snr))
    plot(x(1:size(meanest,2)),meanest(i,:),mycolors{i})
    myleg{end+1} = ['Estimated, snr = ' num2str(snr(i))];
end

if ~isempty(min95est) & ~isempty(max95est)
    plot(x(1:size(min95est,2)),min95est(1,:),'r--');
    plot(x(1:size(min95est,2)),max95est(1,:),'r--');
    myleg{end+1} = ['95% CI for ind. regressions, snr = ' num2str(snr(1))];
end
legend(myleg)


subplot 223; hold on; 
imagesc(msdb,[0 msdblim]), colormap copper
xlabel('Time from trial onset')
ylabel('Signal to Noise Ratio')
%set(gca,'YDir','Reverse')
set(gca,'YTick',1:length(snr))
set(gca,'YTickLabel',snr)
set(gca,'XTick',(1:2:INFO.estlength))
set(gca,'XTickLabel',get(gca,'XTick') .* TR)
title('Mean abs. dev. by time from event onset','FontSize',14)
colorbar('horiz')
    
subplot 224; hold on; 
imagesc(biasmean,bmlim), colormap jet
xlabel('Time from trial onset')
ylabel('Signal to Noise Ratio')
set(gca,'YTick',1:length(snr))
set(gca,'YTickLabel',snr)
set(gca,'XTick',(1:2:INFO.estlength))
set(gca,'XTickLabel',get(gca,'XTick') .* TR)
title('Bias in estimation by time from event onset','FontSize',14)
colorbar('horiz')

return
