%%% plot HRFs for each region as a check %%%%%%%%%%%%%%
for r=1:numreg
    tor_fig(1,size(DATA.fir,2));
    
    for j = 1:size(DATA.fir,2)  % for each event type
        subplot(1,size(DATA.fir,2),j);
        steplot(DATA.fir{r,j});
        title(['Region ' num2str(r) ' Event ' num2str(j)]);
    end
    drawnow
end


%%% plot residuals as a check %%%%%%%%%%%%%%
for s=1:numsub
    for r=1:numreg
        rdat = DATA.resids{s};
        fdat = DATA.fits{s};
        mean_resid(s,:,:) = rdat;
        mean_fit(s,:,:) = fdat;
    end
end

try
    figure;subplot(2,1,1);steplot(mean_resid);title('residuals');
    subplot(2,1,2);steplot(mean_fit);title('fit');
catch, disp('Cannot run steplot');,
end