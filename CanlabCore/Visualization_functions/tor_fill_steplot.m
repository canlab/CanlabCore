function [h,t,dat,d,m1,m2,sterr] = tor_fill_steplot(dat,color,varargin)
% :Usage:
% ::
%
%    [h,t] = tor_fill_steplot(dat,color,[robust flag],[p-thresh],[x vector],[covs no interest])
%
% Plots a mean vector (mean of each column of dat)
% surrounded by a fill with standard err bars
%
% If dat has 3 dimensions, then
% the diff between dat(:,:,1) and dat(:,:,2) is
% used as the difference for computing standard err
% (as in repeated measures)
%
% if behavior is entered as optional argument, removes it before plotting
% lines.  Also returns adjusted output in d, dat
% 
% Optional: robust flag (1/0), robust IRLS
%
% :Examples:
% ::
%
%    tor_fig;
%    tor_fill_steplot(dat,{'b' 'r'},0,.05,secs);

dorobust = 0; pthresh = 0; x = 1:size(dat,2);
t = []; m1 = []; m2 = []; X = [];

if length(varargin) > 0, dorobust = varargin{1};,end
if length(varargin) > 1, pthresh = varargin{2};,end
if length(varargin) > 2, x = varargin{3};,end   % small x, time
if length(varargin) > 3, X = varargin{4};,end   % big X, model matrix of covs

if isempty(x), x = 1:size(dat,2);, end




if length(size(dat)) > 2
    
    % remove covariates of no interest, if any
    % center covs and fit without intercept to preserve mean in data
    if ~isempty(X)
        disp('Adjusting for covariates.');
        X = scale(X,1); % center
        W = X * pinv(X);
        y = squeeze(dat(:,:,1));
        dat(:,:,1) = y - W * y;
    
        y = squeeze(dat(:,:,2));
        dat(:,:,2) = y - W * y;
    end
    
    d = dat(:,:,1) - dat(:,:,2);
    
    if dorobust
        [m1] = robust_mean(dat(:,:,1));
        [m2] = robust_mean(dat(:,:,2));
        [dummy,t,p,sterr] = robust_mean(d);
    else
        
        m1 = nanmean(dat(:,:,1));
        m2 = nanmean(dat(:,:,2));
        md = nanmean(d);
        sterr = ste(d);
        t = md ./ sterr;
        [h,p,ci,stat] = ttest(d);
    end
    
    if ~iscell(color), error('For two groups, color should be cell, e.g.,  {''r'' ''b''})');,end
    
    hold on; 
    h(1) = plot(x,m1,color{1},'LineWidth',2);
    h(2) = plot(x,m2,color{2},'LineWidth',2);
 
    drawnow
    
    fill_around_line(m1,sterr,color{1},x);
    fill_around_line(m2,sterr,color{2},x);

    drawnow
    
    if pthresh
        % significance markers at top of plot.
        yval = max([m1 m2]) + .04 * max([m1 m2]);
        k = 1;
        df = size(d,1) - k;
        %tthr = tinv(1 - pthresh,df);
        tsig = (p < pthresh) .* t;
        
        wh = yval * (tsig > 0);%double((t > tthr)); wh(wh==0) = NaN; wh(wh>0) = yval;
        %plot(x,wh,color{1},'LineWidth',3);
        for i = 1:length(wh)
            if wh(i)
                text(x(i),wh(i),'*','Color',color{1}(1),'FontSize',24);
            end
        end
        
        %wh = yval(find(tsig < 0)); % double((t < -tthr)); wh(wh==0) = NaN; wh(wh>0) = yval;
        %plot(x,wh,color{2},'LineWidth',3);
  
        wh = yval * (tsig < 0);%double((t > tthr)); wh(wh==0) = NaN; wh(wh>0) = yval;
        %plot(x,wh,color{1},'LineWidth',3);
        for i = 1:length(wh)
            if wh(i)
                text(x(i),wh(i),'*','Color',color{2}(1),'FontSize',24);
            end
        end
        
        
        text(x(5),yval+.04 * max([m1 m2]),'Significant','FontSize',16);
    end
    
    drawnow

else
    
    if dorobust
        [m1,t,p,sterr] = robust_mean(dat(:,:,1));
    else
        m1 = nanmean(dat(:,:,1), 1);
        sterr = ste(dat);
    end
    
    if ~iscell(color), tmp = color; color = []; color{1} = tmp; end
    
    hold on; 
    if ischar(color{1})
        h = plot(x,m1,color{1},'LineWidth',2);
    else
        h = plot(x,m1,'o-', 'Color', color{1},'LineWidth',2);
    end
    
    fill_around_line(m1,sterr,color{1},x);
    
end

return
