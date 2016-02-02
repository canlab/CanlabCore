function h = plot_correlation(X,Y,varargin)
% :Usage:
% ::
%
%    handles = plot_correlation(X,Y,varargin)
%
% plots robust or OLS simple or partial correlations
% replaces prplot and plot_correlation_samefig
%
% :Inputs:
%
%   **X:**
%        is matrix of columns of interest plus nuisance
%        default is to plot partial effect of 1st column
%
%   **Y:**
%        is one or more columns of data
%
% :Optional Inputs:
%
%   **'robust':**
%        robust IRLS plot
%
%   **'noprint':**
%        suppress text output
%
%   **'doquad':**
%        quadratic term; not tested, may not work
%
%   **'col':**
%        followed by column of interest
%
%   **'labels':**
%        followed by cell array of text labels for each obs.
%
%   **'colors':**
%        followed by cell array of colors for each column of Y
%
%   **'ylabel':**
%        followed by y-axis label string
%
%   **'xlabel':**
%        followed by x-axis label string
%
%   **'weights':**
%        followed by weights that override any computed ones
%
% :Examples: Plot robust partial corr. 2 of X against col. 17 of Y,
% controlling for other X
% ::
%
%    figure;
%    h = plot_correlation(X,Y(:,17),'col',2,'robust','ylabel','Brain
%                         data','xlabel','Order effect');
%
% Plot Col. 1 of X vs. Y in red squares
% ::
%
%    figure;
%    h = plot_correlation(X,Y(:,17),'robust','colors',{'rs'});
%
%    tor_fig;
%
% ..
%    tor wager, august 2006
% ..

% ---------------------------------------------
% default behaviors
% ---------------------------------------------
dorobust = 0;
doprint = 1;
wh_interest = 1;
doquad = 0;
mylabels = [];
mycol = {'ko' 'rv' 'bs' 'gd' 'y^' 'cv' 'mx'};
ylabelstr = 'Contrast beta value';
xlabelstr = 'Behavioral score';


% ---------------------------------------------
% Optional inputs
% ---------------------------------------------

for i = 1:length(varargin)
    if isstr(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'robust', dorobust = 1;
            case 'noprint', doprint = 0;
            case 'doquad', doquad = 1;
                
            % functional commands
            case 'labels', labels = varargin{i+1};
            case 'colors', mycol = varargin{i+1};
            case 'ylabel', ylabelstr = varargin{i+1}; varargin{i+1} = [];
            case 'xlabel', xlabelstr = varargin{i+1}; varargin{i+1} = [];
            case 'col', wh_interest = varargin{i+1};
            case 'weights', myweights = varargin{i+1};
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% ---------------------------------------------
% intercept
% ---------------------------------------------
wh_intercept = find(all(diff(X) < eps));
if isempty(wh_intercept)
    wh_intercept = size(X,2) + 1;
    X(:,wh_intercept) = 1;
end

ny = size(Y,2);
while length(mycol) < ny, mycol = [mycol mycol]; end

for i = 1:ny
    [x,y,r,p,se,meany,stats] = partialcor(X,Y(:,i),wh_interest,doprint,dorobust);

    % b(1) is param, b(2) is intercept for partial plot
    b = stats.b([wh_interest wh_intercept]);
    
    % weights
    w = stats.w;
    if exist('myweights','var'), w = myweights; end
    
    h{i} = makefigure(x,y,mycol{i},mylabels,doquad,b,w,xlabelstr,ylabelstr);
    
    text(min(X(:,1)),max(y),sprintf('r = %3.2f',r),'FontSize',16,'Color',mycol{i}(1));
    
end


return










function h = makefigure(xvec,yvec,mycol,mylabels,doquad,b,varargin)
% h = makefigure(xvec,yvec,mycol,mylabels,doquad,b,[weights],[xlabel],[ylabel])

ylabelstr = 'Contrast beta value';
xlabelstr = 'Behavioral score';

hold on; grid on; set(gca,'FontSize',18)

if length(varargin) > 0
    % ROBUST (or we just have weights)
    w = varargin{1};
    for i = 1:length(xvec)
        h(i) = plot(xvec(i),yvec(i),mycol,'LineWidth',.5,'MarkerSize',8, ...
            'MarkerFaceColor',mycol(1));
        set(h(i),'MarkerFaceColor',[repmat(1-w(i),1,3)] )
    end

    % set axis
    xlims = [min(xvec) max(xvec)];
    ylims = [min(yvec) max(yvec)];
    xrange = (xlims(2) - xlims(1)) *.1; % % margin
    yrange = (ylims(2) - ylims(1)) *.1;
    xlims = xlims + [-xrange xrange];
    ylims = ylims + [-yrange yrange];
    
    set(gca,'Xlim',xlims,'YLim',ylims);

    % plot regression line by hand
    xl = get(gca,'Xlim'); yl = get(gca,'Ylim');
    x = xl(1)-5:xl(2)+5;
    h2 = plot(x,b(1) * x + b(2),'-','Color',mycol(1),'LineWidth',2);
    h = [h h2];
    set(gca,'XLim',xl,'YLim',yl)


else
    % not robust

    h = plot(xvec,yvec,mycol,'LineWidth',3,'MarkerSize',6,'MarkerFaceColor',mycol(1));

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

end


drawnow

if length(varargin) > 1
    xlabelstr = varargin{2};
end

if length(varargin) > 2
    ylabelstr = varargin{3};
end

if ~isempty(mylabels)
    for j = 1:length(xvec)
        text(xvec(j),yvec(j),mylabels{j},'FontWeight','b')
    end
end

ylabel(ylabelstr)
xlabel(xlabelstr)

return
