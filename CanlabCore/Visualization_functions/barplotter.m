function h = barplotter(data, varargin)
% :Usage:
% ::
%
%    h = barplotter(data)
%
% Creates a barplot out of the columns or elements (if data is a cell
% array containing vectors) of data.
% ::
%
%    h =barplotter(..., 'groups', grouping)
%
% grouping must have the same number of columns or elements (if data is a
% cell array) as data, and must consist of integers beginning with 1 and
% ending with the total number of groups. Data belonging to the same group
% will be 'clustered' together in the plot. Ommitting groups (eg, [1 1 3 3])
% will create extra spacing between groups.
% ::
%
%    h =barplotter(..., 'std')
%
% overrides the default behavior and plots standard deviation bars instead
% of standard error bars
% ::
%
%    h =barplotter(..., 'CI', alpha)
%
% this will override the default behavior (as well as 'std') and plot
% 1-alpha confidence intervals.
% ::
%
%    h =barplotter(..., 'labels', labels)
%
% labels must be a cell array of strings with the same number of elements
% as data has columns or elements. The x-axis will be labeled with these.
% ::
%
%    h =barplotter(..., 'label_groups')
%
% Applies labels to groups instead to individual bars
% ::
%
%    h =barplotter(..., 'legend', names)
%
% Plots a legend in the figure, with labels corresponding to the elements
% of the cell vector of strings, names.
% ::
%
%    h =barplotter(..., 'plegend', names, p)
%
% As 'legend', above, but p is a vector of the barplots to include in the
% legend (e.g. if you plotted 6 bars and p = [1 3], only the first and third
% bar would be included in the legend.
% ::
%
%    h =barplotter(..., 'PlotLineHor', value)
%
% Plots a horizontal line at the Y value indicated.
% 
% h = barplotter(..., 'LinePlot', colors, styles, markers)
% 
% Changes the graph from a bar graph to a line graph. Setting groups will
% cause one line to be drawn for each group.
% 
% h = barplotter(..., 'ErrorWidth', errorwidth)
% 
% Sets the line weight (in points) for error bars
% ::
%
%    h =barplotter(..., 'PropertyName', PropertyValue)
%
% Properties correspond to various Matlab figure properties, as
% appropriate. Currently supported properties (more to be added) are:
%
%   - 'Title'
%   - 'XLabel'
%   - 'YLabel'
%   - 'YLim'
%   - 'XLim'
%   - 'YTick'
%   - 'YTickLabel'
%   - 'YMinorTick'
%   - 'FontSize'
%   - 'xFontSize' %for xlabel
%   - 'yFontSize' %for ylabel
%   - 'tFontSize' %for title
%   - 'FaceColor' %note that you may specify a matrix of 3-element RGB vectors,
%               rather than a single vector. Barplotter will then cycle
%               the rows of the matrix until all bars have been drawn.
%   - 'Colormap'
%   - 'GridLineStyle' % - | - -| {:} | -. | none
%   - 'TickDir' % in or out
%   - 'MarkerSize'
%   - 'LineWidth'
%
% ..
%    Note that this function has no native error handling. If you're
%    encountering inexplicable errors it's likely because you haven't passed
%    in the correct inputs.
%    02/10/07 Jared Van Snellenberg
% ..


if ~iscell(data)
    X = data;
    clear data
    for k = 1:size(X, 2)
        data(k) = {X(:,k)};
    end
    clear X
else
    data = data(:);
end

SD = 0;
CI = 0;
fc = 'flat';
grouping = ones(size(data))';
colormapfun = @hsv;
xlabels = [];
Title = [];
Xlabel = [];
Ylabel = [];
fontsize = 18;
xfontsize = 20;
yfontsize = 20;
tfontsize = 20;
lgroups = 0;
leg = 0;
grd='none';
drawline=0;
lineplot=0;
yminortick='off';
tickdir = 'out';
markersize = 12;
linewidth = 2;
errorwidth = 1;

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch(varargin{k})
            case 'groups'
                grouping = varargin{k+1};
            case 'std'
                SD = 1;
            case 'CI'
                CI = 1;
                alpha = varargin{k+1};
            case 'labels'
                xlabels = varargin{k+1};
            case 'label_groups'
                lgroups = 1;
            case 'YLim'
                yl = varargin{k+1};
            case 'XLim'
                xl=varargin{k+1};
            case 'YTick'
                ytick = varargin{k+1};
            case 'YTickLabel'
                yticklabels=varargin{k+1};
            case 'YMinorTick'
                yminortick='on';
            case 'FaceColor'
                fc = varargin{k+1};
            case 'Colormap'
                colormapfun = str2func(varargin{k+1});
            case 'Title'
                Title = varargin{k+1};
            case 'XLabel'
                Xlabel = varargin{k+1};
            case 'YLabel'
                Ylabel = varargin{k+1};
            case 'FontSize'
                fontsize = varargin{k+1};
            case 'xFontSize'
                xfontsize = varargin{k+1};
            case 'yFontSize'
                yfontsize = varargin{k+1};
            case 'tFontSize'
                tfontsize = varargin{k+1};
            case 'legend'
                leg = 1;
                names = varargin{k+1};
                p = 1:length(data);
            case 'plegend'
                leg = 1;
                names = varargin{k+1};
                p = varargin{k+2};
            case 'GridLineStyle'
                grd=varargin{k+1};
            case 'PlotLineHor'
                drawline=1;
                yLine=varargin{k+1};
            case 'LinePlot'
                lineplot = 1;
                colors = varargin{k+1};
                styles = varargin{k+2};
                markers = varargin{k+3};
            case 'TickDir'
                tickdir = varargin{k+1};
            case 'MarkerSize'
                markersize = varargin{k+1};
            case 'LineWidth'
                linewidth = varargin{k+1};
            case 'ErrorWidth'
                errorwidth = varargin{k+1};
        end
    end
end

if isnumeric(fc)
    while size(fc, 1) ~= length(data)
        fc = [fc;fc];
        if size(fc, 1)>length(data)
            fc(length(data)+1:end,:) = [];
        end
    end
    cm = fc;
    fc = 'flat';
else
    cm = colormapfun(length(data));
end
if strcmp(fc, 'none')
    for k = 1:length(data)
        C{k} = 'none';
    end
else
    for k = 1:length(data)
        C{k} = cm(k,:);
    end
end


for k = 1:length(data)
    means(k) = mean(data{k});
    if CI
        err(:,k) = [means(k)+tinv(1-alpha/2, length(data{k}-1))*std(data{k})/sqrt(length(data{k}))...
            means(k)-tinv(1-alpha/2, length(data{k}-1))*std(data{k})/sqrt(length(data{k}))];
    elseif SD
        err(:,k) = [means(k)+std(data{k}); means(k)-std(data{k})];
    else
        err(:,k) = [means(k)+std(data{k})/sqrt(length(data{k})); means(k)-std(data{k})/sqrt(length(data{k}))];
    end
end

if ~exist('xl','var')
    if ~lineplot
        xl = [0 max(grouping)*0.5+length(data)+0.5];
    else
        xl = [0 (length(data) / max(grouping) + 1)];
end

if ~exist('yl', 'var')
    yl = max(max(err));
    if yl>0
        count = 0;
        if yl>1
            while floor(yl)
                yl = yl/10;
                count = count+1;
            end
            yl = yl*10;
            count = count-1;
        else
            while ~floor(yl)
                yl = yl*10;
                count = count-1;
            end
        end
        if yl>0.9
            yl = ceil(yl)*10^count;
        else
            yl = ceil(yl*1.1)*10^count;
        end
        yl = [min([0 min(min(err))+min(min(err))/10]) yl];
    else
        yl = min(min(err));
        count = 0;
        if yl < -1
            while ceil(yl)
                yl = yl/10;
                count = count+1;
            end
            yl = yl*10;
            count = count-1;
        else
            while ~ceil(yl)
                yl = yl*10;
                count = count-1;
            end
        end
        if yl < -0.9
            yl = floor(yl)*10^count;
        else
            yl = floor(yl*1.1)*10^count;
        end
        yl = [yl 0];
    end
end

figure;
h = axes('GridLineStyle',grd,'YGrid','on','YMinorTick',yminortick);
plotcount = 0;
Xtick = [];
set(h, 'XLim', xl, 'YLim', yl)
if drawline
    hold on
    line(xl,[yLine yLine],'Color','k','LineStyle','--','LineWidth',2);
    hold off
end

if lineplot
    for k = 1:max(grouping)
        plotcount = plotcount+1;
        grp = find(grouping == k);
        if size(grp,1)>1
            grp=grp';
        end
        hold on
        A(plotcount) = line(1:length(grp),means(grp),'Color',colors{k},'LineStyle',styles{k},'Marker',markers{k},'MarkerSize',markersize,'LineWidth',linewidth);
        count=1;
        for j = grp
            line([count count-0.1 count-0.1; count count+0.1 count+0.1],[err(1, j) err(1, j) err(2, j); err(2, j) err(1, j) err(2, j)],'Color',colors{k},'LineStyle',styles{1},'LineWidth',errorwidth);
            count = count + 1;
        end
        hold off
    end
    Xtick = 1:sum(grouping == 1);
else
    for k = 1:max(grouping)
        grp = find(grouping == k);
        if size(grp,1)>1
            grp=grp';
        end
        for j = grp
            plotcount = plotcount+1;
            hold on
            A(plotcount) = area(h, [grouping(j)*0.5+plotcount-0.5-0.4 grouping(j)*0.5+plotcount-0.5+0.4], [means(j) means(j)], 'FaceColor', C{j});
            line(...
                [grouping(j)*0.5+plotcount-0.5 grouping(j)*0.5+plotcount-0.5-0.2 grouping(j)*0.5+plotcount-0.5-0.2;...
                grouping(j)*0.5+plotcount-0.5 grouping(j)*0.5+plotcount-0.5+0.2 grouping(j)*0.5+plotcount-0.5+0.2], ...
                [err(1, j) err(1, j) err(2, j); err(2, j) err(1, j) err(2, j)], 'Color', 'k','LineWidth',errorwidth);
            hold off
            if lgroups
                if j == grp(end)
                    Xtick(end+1) = grouping(j)*0.5+plotcount-length(grp)/2;
                end
            else
                Xtick(end+1) = grouping(j)*0.5+plotcount-0.5;
            end
        end
    end
end
if leg
    legend(A(p), names{:});
end
set(h, 'XTick', Xtick, 'XTickLabel', xlabels)
if exist('ytick', 'var')
    set(h, 'YTick', ytick)
end
if exist('yticklabels','var')
    set(h,'YTickLabel',yticklabels)
end
xlabel(h, Xlabel, 'FontSize', xfontsize)
ylabel(h, Ylabel, 'FontSize', yfontsize)
title(h, Title, 'FontSize', tfontsize)
set(h,'FontSize', fontsize)
set(h,'TickDir',tickdir)
end
