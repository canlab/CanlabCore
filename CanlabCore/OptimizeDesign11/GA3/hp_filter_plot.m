function [allxc,allip,allgp,allzip,allzgp] = hp_filter_plot(HP,acfscale,X,TR,varargin)
%function [allxc,allip,allgp,allzip,allzgp] = hp_filter_plot(HP,acfscale,X,TR,[tsim])
%
% Tor Wager, 2/25/04
% Saves and plots power estimates as a function of high-pass filtering
% length and autocorrelation parameters
%
% HP is a vector of lengths (e.g., 20:10:300)
% acfscale is a vector of autocorrelation scalings (canonical * acfscale)
% e.g., [2 1 .5 .3 .2 .1]
%
% tsim = optional argument, 1 indicating that you want to run actual
% regressions with simulated data to compare to predictions from algebra.
%
% example: sac is onset times for saccade task, TR = 2
% TR = 2; hrf2 = spm_hrf(TR) ./ max(spm_hrf(TR)); 
%HP = 20:10:200; acfscale = [2 1 .5 .3 .2 .1];
% [Xactual,delta,delta_hires,hrf] = onsets2delta(sac,TR);
% dhr = getRandom(delta_hires);
% [tmp,d] = downsample_delta(dhr,16*TR); X=getPredictors(d,hrf2);
% X(:,end+1) = 1;
% [allxc,allip,allgp] = hp_filter_plot(HP,acfscale,X,TR);

warning off
time1 = clock;

if length(varargin) > 0, tsim = varargin{1};, f1=figure('Color','w');, else, tsim = 0;, end

con = [1 -1];   % contrast matrix
n_in_group = 15;

f0 = figure('Color','w'); colors = {'r' 'g' 'b' 'k' 'y' 'c' 'm'};

if length(acfscale) > length(colors), error('Choose no more than 7 acfscale values'),end

aind = 1;


for acf = acfscale
    
    ind = 1; clear ip, clear gp
    [xc2,Vi] = canonical_autocorrelation(TR,size(X,1),acf);
    allxc{aind} = xc2;
    
    subplot(1,3,1); hold on ; set(gca,'FontSize',16)
    plot(allxc{aind},colors{aind}), drawnow

    for HPlength = HP

        [S,Vi] = getSmoothing(HPlength,0,TR,size(X,1),xc2);
        
        [t1,t2] = xpower(X,con,[],n_in_group,[],Vi,S);
        
        ip(ind) = t1(1);
        gp(ind) = t2(1);
    
        % optional tsim
        if tsim
            [izval,gzval,OUT] = xzpower(X,con,[],n_in_group,[],Vi,S);
            zip(:,ind) = izval;
            zgp(:,ind) = gzval;
        end
            
        ind = ind + 1;
    
    end
    
    h(1) = subplot(1,3,2); hold on; set(gca,'FontSize',16),plot(HP,ip,colors{aind}); 
    xlabel('High-pass filter length'), ylabel('Power (Z)'),;title('Individual power')
    h(2) = subplot(1,3,3); hold on; set(gca,'FontSize',16),plot(HP,gp,colors{aind}); 
    xlabel('High-pass filter length'), ylabel('Power (Z)'),;title('Group power, n = 15')

    drawnow
    
    if tsim
        figure(f1); 
        hi(1) = subplot(2,2,1); hold on; set(gca,'FontSize',16),plot(HP,ip,colors{aind}); 
        xlabel('High-pass filter length'), ylabel('Power (Z)'),;title('Individual power')
        hi(2) = subplot(2,2,2); hold on; set(gca,'FontSize',16),plot(HP,zip,colors{aind}); 
        xlabel('High-pass filter length'), ylabel('Power (Z)'),;title('Individual power (full sim)')
        hi(3) = subplot(2,2,3); hold on; set(gca,'FontSize',16),plot(HP,gp,colors{aind}); 
        xlabel('High-pass filter length'), ylabel('Power (Z)'),;title('Group power, n = 15')
        hi(4) = subplot(2,2,4); hold on; set(gca,'FontSize',16),plot(HP,zgp,colors{aind}); 
        xlabel('High-pass filter length'), ylabel('Power (Z)'),;title('Group power (full sim)')
        drawnow
        figure(f0)
        
        allzip(aind,:) = zip;
        allzgp(aind,:) = zgp;
    end
    
    allip(aind,:) = ip;
    allgp(aind,:) = gp;
    
    aind = aind + 1;
    
end

warning on

equalize_axes(h)
subplot(1,3,1)
for i  = 1:length(acfscale), str{i} = sprintf('%3.2f',acfscale(i));, end

legend(str)

    if tsim
        figure(f1); equalize_axes(hi);
    end

    t2 = etime(clock,time1);
    disp(['Finished in ' num2str(t2) ' s'])
    
return