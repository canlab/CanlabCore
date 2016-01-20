function [rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,ALLINFO,hrf,snr,TR] = ideal_deconv6(conditions,mspec,ttype)
% Tests deconvolution matrix directly against idealized data
% you put in the exact temporal sequence to be deconvolved,
% in the form of the DX matrix.
%
% :Usage:
% ::
%
%     [rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,ALLINFO,hrf,snr,TR] = ideal_deconv6(conditions,mspec,ttype)
%
% :Inputs:
%
%   **DX:**
%        deconvolution matrix
%
%   **tp:**
%        time points estimated for each condition in DX
%
%   **TR:**
%        repetition time of scan
%
%   **ttype:**
%        trial types to test (out of 1:n different conditions in DX)
%        recommended for time saving to use ttype = a single number only
%
% This function is like ideal_deconv5, but tests variability across designs as well.
%
% ..
%    Tor Wager, 4/19/02
% ..

hrf = spm_hrf(mspec.TR);
hrf = hrf ./ max(hrf);

% hrfin contains cell array of hrfs for each condition
[DX] = construct_model(mspec,conditions,[]);
[delta, wd, wb] = DX_find_delta(DX);
for j = 1:length(wd)
    hrfin{j} = hrf;
    
    % when epoch or block regressor - now only works for one-part designs
    if length(conditions) >= j
        if isfield(conditions(j),'stimlength')
            if conditions(j).stimlength > 1
                hrfin{j} = conv(hrf,ones(conditions(j).stimlength,1));
                hrfin{j} = hrfin{j} ./ max(hrfin{j});
            end
        end
    end
    
end

% for custom hrf for one trial type
%hrf2 = sum([[zeros(5,1);hrf] [hrf;zeros(5,1)]],2);
%hrfin{2} = hrf2;


ndesigns = 2;
nnoise = 100;
snr = 1;
designs = 1:ndesigns;

for j = ttype
    
    clear rmsd, clear msdstd, clear msdb, clear biasmean, clear meanest, clear min95est, clear max95est
    
    for i = 1:ndesigns   % for this many designs
    
    % rmsd      sqrt of mean squared dev. of estimates from true response
    % msdstd    single obs. 95% confidence interval for error (2-tailed)
    % msdb      abs dev. from ideal response for each time point
    % biasmean  mean bias (above or below true response) for each time pt.
    % meanest   mean estimate of the response (estimated hrf)
    % min95est  95% of individual regressions fall within this window...
    % max95est  between min95est and max95est

        [DX] = construct_model(mspec,conditions,[]);
        
        [rmsd(i),msdstd(i,:),msdb(i,:),biasmean(i,:),meanest(i,:),min95est(i,:),max95est(i,:),INFO] ...
        = dx_estimate_dev(DX,snr,hrfin,j,nnoise);
    
    end

    
    
    % plotting
    
    
    
    plot_ideal_deconv5(rmsd,msdstd,msdb,biasmean,meanest,min95est,max95est,INFO,hrfin{j},designs,mspec.TR);

    % re-plot first graph to fit this format
    subplot 221
    cla
    plot([(1:ndesigns); (1:ndesigns)], msdstd','k-','LineWidth',2)
    hold on;
    plot(rmsd,'k.','LineWidth',2)
    xlabel('Design Realization')

    subplot 222
    cla
    plot(repmat((1:size(meanest,2))' .* mspec.TR,1,size(meanest,1)), meanest')
    h = plot((1:length(hrf)) * mspec.TR, hrf,'k','LineWidth',2);
    if ~isempty(INFO.autocorrelation)
        h(2) = plot((1:length(INFO.autocorrelation))./(length(INFO.autocorrelation) / length(hrf)),INFO.autocorrelation,'k--','LineWidth',2);
        myleg = {'HRF','ACF'};
    else
        myleg = {'HRF - white noise simulated'};
    end
    legend(h,myleg)

    INFO.TR = mspec.TR;
    ALLINFO(j) = INFO;

    a = get(gcf,'Children');
    axes(a(3))
    ylabel('Design Realization')
    set(gca,'YTick',[1 ndesigns])
    set(gca,'YTickLabel',[1 ndesigns])
    
    %subplot 223
    %hold on
    %ylabel('Design Realization')
    %set(gca,'YTickLabel',[1 ndesigns])
    
    axes(a(5))
    ylabel('Design Realization')
    set(gca,'YTick',[1 ndesigns])
    set(gca,'YTickLabel',[1 ndesigns])

    %subplot 224
    %hold on
    %ylabel('Design Realization')
    %set(gca,'YTickLabel',[1 ndesigns])
   
    figure;set(gcf,'Color','w')
    a = mean((biasmean),2);
hist(a,length(a) / 4)
title(['Frequency histogram of mean estimation bias across ' num2str(nnoise) ' noise realizations per model'])
end

return
