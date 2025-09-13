function [r,str,sig,ploth,rp] = plot_correlation_samefig(xvec,yvec,varargin)
% :Usage:
% ::
%
%    [r,infostring,sig,h] = plot_correlation_samefig(xvec,yvec,[textlabs],[color],[doquad],[dorobust])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **xvec:**
%        x vector
%
%   **yvec:**
%        y vector
%
%   **varargin**
%        string of text labels
%
% :Optional Inputs:
%   **doquad:**
%        flag for quadratic correlations as well!
%
%   **dorobust:**
%        remove n outliers from data, using Min Cov Determinant (MCD)
%        Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
%        Journal of the American Statistical Association, Vol. 79, pp. 871-88
%        outliers calculated using IRLS (robustfit.m) do not work well.
%        you enter n
%
%   empty variable arguments are OK, defaults will be used
%
% :Examples:
% ::
%
%    figure; [r, infos] = plot_correlation_samefig(x, y, [], 'ko', 0, 1);
%
%    % for text labels only, try:
%    plot_correlation(beh1,mri1,highlow,'w.');


robopt = 'IRLS';    %'IRLS' or 'MCD'

xvec = double(xvec);
yvec = double(yvec);

r = [];, str=[];sig=[]; hh=[];, p=[];

if length(varargin) > 0, mylabels = varargin{1};  else  mylabels = [];  end
if length(varargin) > 1, mycol = varargin{2};  else  mycol = 'ko'; end
if length(varargin) > 2, doquad = varargin{3};  else  doquad = 0;  end
if length(varargin) > 3, dorobust = varargin{4};  else  dorobust = 0;  end
if isempty(mycol), mycol = 'ko'; end

if size(xvec,2) > size(xvec,1) && size(xvec,1)==1, xvec = xvec';  end
if size(yvec,2) > size(yvec,1) && size(yvec,1)==1, yvec = yvec';  end


wh = find(isnan(xvec) | isnan(yvec));
xvec(wh,:) = []; yvec(wh) = [];


X = [xvec ones(length(xvec),1)];
y = (yvec);

% robust regression - remove outliers
if dorobust && strcmp(robopt,'MCD')
    
    %if doquad
    %tmp = xvec .^2; tmp = tmp-mean(tmp);
    %[res]=fastmcd_noplot([tmp X(:,1) y]);
    %else
    [res]=fastmcd_noplot([X(:,1) y]);
    %end
    
    % remove n most extreme outliers and recompute correlation
    
    wh = res.flag==0; nout = sum(res.flag==0);
    y(wh) = []; X(wh,:) = []; xvec(wh) = []; yvec(wh) = [];
    
    %tmp = [(1:length(y))' stats.w]; tmp=sortrows(tmp,2); wh=tmp(1:dorobust,1);
    
    %keyboard
    %for i = 1:length(y),lab{i}=num2str(stats.w(i)); end
    %    makefigure(xvec,yvec,'wo',lab,doquad)
    %    hold on; plot(xvec(wh),yvec(wh),'ro')
    %end
    ploth = makefigure(X(:,1),y,mycol,mylabels,doquad,b);
    
    
elseif dorobust && strcmp(robopt,'IRLS')
    
    % to get b and stats in original scale
    [b,stats]=robustfit(X,y,'bisquare',[],'off');
    rZ = stats.t(1); rp = stats.p(1);
    if rp < .05, sig=1; else sig=0; end
    
    % correlation (r) is NOT the beta of standardized variables, as
    % with OLS, because the standardized beta doesn't account for the
    % change in covariance due to the weights.
    [r] = weighted_corrcoef([X(:,1) y],stats.w);
    r = r(1,2);
    
    ploth = makefigure(X(:,1),y,mycol,mylabels,doquad,b,stats.w);
    
else
    % ols
    b = X \ y;
    ploth = makefigure(X(:,1),y,mycol,mylabels,doquad,b);
end




% text and stats




if dorobust & strcmp(robopt,'IRLS')
    
    str = sprintf('Sig. of B0: u=%3.2f, t=%3.2f, p=%3.4f\n C: u=%3.2f, t=%3.2f, p=%3.4f   R: r=%3.2f, t=%3.2f, p=%3.4f', ...
        mean(y),stats.t(end),stats.p(end),b(1), stats.t(1), stats.p(1), r, rZ, rp);
    
    %text(min(X(:,1)),max(y),sprintf('r = %3.2f',r),'FontSize',16)
    
    
else
    % not robust IRLS
    
    % in regression, test significance of intercept parameter
    y2 = (y - X*b);   % subtract b1, the regression fit
    y2 = y2 + mean(y);  % add intercept back in
    % Old, not necessary to do it this way...
    
    try
        [h,b0p,ci,b0stats] = ttest(y2);
    catch
        disp('No ttest.m: No stats toolbox?')
    end
    
    
    % separate correlation and t-test
    try
        [h,p,ci,stats] = ttest(y);
    catch
        disp('No ttest.m: No stats toolbox?')
    end
    
    r = corrcoef(y,X(:,1)); r= r(1,2);
    pr = partialcorr(y,X(:,1),X(:,2:end));

    try
        [rci,sig,rZ,rp] = r2z(r,length(y),.05);
    catch
        disp('Error in or missing r2z.m: No stats toolbox?')
    end
    
    ht=text(min(X(:,1)),max(y),sprintf('r = %3.2f',pr),'FontSize',20);
    ploth = [ploth ht];
    
    try
        str = sprintf('Sig. of B0: u=%3.2f, t=%3.2f, p=%3.4f\n C: u=%3.2f, t=%3.2f, p=%3.4f   R: r=%3.2f, Z=%3.2f, p=%3.4f', ...
            mean(y2),b0stats.tstat,b0p,mean(y), stats.tstat, p, r, rZ, rp);
    catch
        str = ['Missing stats'];
        sig = NaN;
    end
    if dorobust && strcmp(robopt,'MCD'), str=[str sprintf(' N_o_u_t=%3.0f',nout)]; end
    
end % Robust or not

title(str,'FontSize',12)


if doquad
    tmp = xvec .^ 2; tmp = tmp - mean(tmp);
    X = [tmp xvec ones(length(xvec),1)];
    X(isnan(xvec),:) = [];
    y = (yvec);
    y(isnan(xvec)) = [];
    
    [b,bint,dummy,rint,stats] = regress(y,X,.05);
    
    r0 = r;
    r = sqrt(stats(1));
    
    mysig = sum(sign(bint),2);
    for i = 1:size(b,1),if mysig(i),sigstr{i}='*'; else  sigstr{i}=''; end, end
    
    str = sprintf('y = %3.2fX^2%s + %3.2fX%s + %3.2f%s, r_0=%3.2f, Mult. r=%3.2f, Omni F=%3.2f, p=%3.2f', ...
        b(1),sigstr{1},b(2),sigstr{2},b(3),sigstr{3},r0, r, stats(2),stats(3));
    
    sig = stats(3) < .05;
    if dorobust & strcmp(robopt,'MCD'), str=[str sprintf(' N_o_u_t=%3.0f',nout)]; end
    
    title(str,'FontSize',12)
    refcurve(b)
    
end

end % function





function h = makefigure(xvec,yvec,mycol,mylabels,doquad,b,varargin)

hold on; grid off; set(gca,'FontSize',18)

if length(varargin) > 0
    % ROBUST
    w = varargin{1};
    for i = 1:length(xvec)
%         h(i) = plot(xvec(i),yvec(i),[mycol(1) 'o'],'LineWidth',.5,'MarkerSize',8, ...
%             'MarkerFaceColor',mycol(1));
          h(i) = scatter(xvec(i),yvec(i),50,mycol(1),'o','filled');
%         set(h(i),'MarkerFaceColor',[repmat(1-w(i),1,3)] )
%         set(h(i), 'Color', [.2 .2 .2], 'LineWidth',1)
          
    end
    
    % set axis
    set_axis(xvec, yvec);
    
    % plot regline by hand
    xl = get(gca,'Xlim'); yl = get(gca,'Ylim');
    x = xl(1)-5:xl(2)+5;
    plot(x,b(1) * x + b(2),'k-','LineWidth',.5)
    set(gca,'XLim',xl,'YLim',yl)
    
    
else
    % not robust
    
%     h = plot(xvec,yvec,mycol,'LineWidth',3,'MarkerSize',6,'MarkerFaceColor',mycol(1));
    h = scatter(xvec,yvec,50,mycol,'o','filled');
%     set(h, 'Color', [.2 .2 .2], 'LineWidth',1)
    
    if doquad
        %refcurve(doquad)
    else
        try
            refline
        catch
            tmp = get(gca,'XLim');
            x = min(xvec)-std(xvec):min(std(xvec),.01):max(xvec)+std(xvec);
            plot(x,b(1) * x + b(2),[mycol(1) '-'],'LineWidth',.5)
            set(gca,'XLim',tmp);
            
        end
    end
    
    set_axis(xvec, yvec);

    
end


drawnow



if ~isempty(mylabels)
    for j = 1:length(xvec)
        ht = text(xvec(j),yvec(j),mylabels{j},'FontWeight','b');
        h = [h ht];
    end
end

ylabel('Contrast beta value')
xlabel('Behavioral score')

end %function




function set_axis(xvec, yvec)

% set axis
xlims = [min(xvec) max(xvec)];
ylims = [min(yvec) max(yvec)];
xadj = range(xlims) .* [-.1 .1];
if xadj(1) == 0, xadj(1) = -.1;  end
if xadj(2) == 0, xadj(2) = .1;  end
xlims = xlims + xadj;

yadj = range(xlims) .* [-.1 .1];
if yadj(1) == 0, yadj(1) = -.1;  end
if yadj(2) == 0, yadj(2) = .1;  end
xlims = xlims + yadj;
ylims = ylims + yadj;
set(gca,'Xlim',xlims,'YLim',ylims);

end
