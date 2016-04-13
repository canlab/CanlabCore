function hewma2_plot(Zpop,mu,cp,cp_ind,ooc_indices,maxtime,tthresh,tt,tm,q,Wm,Zm,sterr,varargin)
% :Usage:
% ::
%
%     hewma2_plot(Zpop,mu,cp,cp_ind,ooc_indices,maxtime,tthresh,tt,tm,q,Wm,Zm,sterr,varargin)
%
% see hewma2.m

cp = round(cp);

dogroupplot = 0;
tagname  = 'hewma2display';

if length(varargin) > 0
    dogroupplot = 1;
    tagname  = 'hewma2diffdisplay';
end


f1 = findobj('Tag', tagname);
if isempty(f1)  
    f1 = tor_fig(2,2); 
    set(f1, 'Tag', tagname);
else
    figure(f1)
    for i = 1:4, subplot(2, 2, i); cla; end
end

ZZ = Zpop-mu;

% ----------------------------------------------------------
% *
% * PLOT 1: GROUP AVERAGE (OR DIFFERENCE)
% *
% ----------------------------------------------------------
subplot(2, 2, 1);

hh(1) = plot(ZZ,'k','LineWidth',2);

% baseline box
% ------------------------------
    yl = get(gca,'YLim');
    cl = tthresh .* mean(sterr);
    yl(1) = min([yl(1)-1 -cl-1]); yl(2) = max([yl(2)+1 cl+1]);
    
    fill([1 tt tt 1],[yl(1) yl(1) yl(2) yl(2)],[.85 .85 .85]);
    text(5,min(ZZ),'Baseline','FontSize',16);
    
% group Z stat
% ------------------------------   
hh(1) = plot(ZZ,'k','LineWidth',2);
hold on

 % change points
% ------------------------------
    if ~isnan(cp)
        plot(cp,ZZ(cp),'bo','MarkerSize',8,'MarkerFaceColor','g');
        plot([cp cp],get(gca,'YLim'),'g','LineWidth',2);
    end
    
    % individual change points
    if ~isempty(cp_ind)
        if ~any(isnan(cp_ind))
            plot(cp_ind,ZZ(cp_ind),'go','MarkerSize',6);
        end
    end
 
% OOC points in red
% ------------------------------
    % all ooc points
    ZZtmp = NaN * zeros(size(ZZ)); ZZtmp(ooc_indices) = ZZ(ooc_indices);
    plot(ZZtmp,'r','LineWidth',2);
  
% max t value in blue
% ------------------------------

    if ~isnan(maxtime)
        plot(maxtime,ZZ(maxtime),'go','MarkerSize',8,'MarkerFaceColor','b');
    end
   
 % standard error
% ------------------------------
    fill_around_line(ZZ,sterr,'k');  
    
    
    title('Group activation'); 
    if length(varargin) > 0, title('Difference between groups'); end 
    
 % control limits
% ------------------------------   
    cl = tthresh .* mean(sterr);
    hh(2) = plot(get(gca,'XLim'),[cl cl],'k--');
    plot(get(gca,'XLim'),[-cl -cl],'k--');
    legend(hh,'Group mean','Control limit (approximate)');
    
    
    
% ----------------------------------------------------------
% *
% * PLOT 2: NULL HYPOTHESIS DISTRIBUTION
% *
% ----------------------------------------------------------   
    
    subplot(2,2,2);
    hist(q,50); hold on; plot([abs(tm) abs(tm)],get(gca,'YLim'),'k-','LineWidth',2);
    title('Null-hypothesis max t-value and observed max t-value');

% ----------------------------------------------------------
% *
% * PLOT 3: SUBJECT WEIGHTS
% *
% ----------------------------------------------------------

    subplot(2,2,3);
    n = 1:length(Wm);
    plot([n;n],[Wm';Wm'],'o','LineWidth',2); 
    title('Case weights');
 
% ----------------------------------------------------------
% *
% * PLOT 4: INDIVIDUAL SUBJECTS (OR GROUPS) 
% *
% ----------------------------------------------------------

    subplot(2,2,4);
    
    % group contrast, or individuals
    if dogroupplot
        gfits = varargin{1};
        plot(gfits'); 
        title('Groups');
        legend({'Low' 'High'});
        yval = max(gfits(:)) + .01*max(gfits(:));
        %yval = repmat(yval,1,length(ZZtmp));
        ZZtmp(~isnan(ZZtmp)) = yval;
        hold on;
        
        plot(ZZtmp,'r','LineWidth',3);
        
        if ~isnan(cp)
            plot([cp cp],get(gca,'YLim'),'g','LineWidth',2);
        end
        
        % baseline box
        yl = get(gca,'YLim');
        fill([1 tt tt 1],[yl(1) yl(1) yl(2) yl(2)],[.85 .85 .85]);
        text(5,min(ZZ),'Baseline','FontSize',16);
    
    else
        plot(Zm'); 
        title('Individual cases');
        drawnow
    end
    
    
    return
