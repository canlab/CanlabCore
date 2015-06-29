function [h, axh]=barplot_colored(data,varargin)

% haven't had time to document this yet
% this is a good function though
% within-subject error bars now added; use 'within'
%
% [bar_handles, axis_handle]=barplot_colored(data,varargin)
% 
% Input arguments: Optional 
% 'within' : Do within-subject STE bars, average obs x variable interaction
%            Loftus and Masson 1994 style.
%
% Strings, followed by values for each:
%       < Color control >
% 'colormap' : followed by colormap name to use
% 'colors' : followed by cell array of colors per bar; supercedes colormap
%
%       < Display items >
% 'fontsize'
% 'title'
% 'XTickLabels'
% 'ylabel'
% 'xlabel'
%
%       < bar locations >
% 'x' : followed by x values for bars (locations)

%       
% NOTE: For this function, keywords must be even-numbered argument entries,
% e.g., arg 2, 4, 6.  Odd argument entries are values.
% For example: This works, and you need the extra empty arg after 'within'
% [h1, s1] = barplot_colored(pexp1, 'within', ' ', 'title', 'Pattern expression', 'XTickLabels', dat.Y_names, 'x', 1:nterms);
%
% You can assign arbitrary colors to bars by setting the colormap:
% [h, s] = barplot_colored([corr_temp corr_rep]);
% cm = [1 .5 0; .5 0 1];
% colormap(cm)
%
% Example: A grouped barplot
% ---------------------------------------------------
% dat = rand(20, 4);
% create_figure('bars');
% [h1, s1] = barplot_colored(dat, 'x', [1 2 4 5]);
% % set(h2, 'BarWidth', .9)
%
% Change colormap:
%[h1, s1] = barplot_colored(dat, 'x', [1 2 4 5], 'colormap', 'summer');
%
% Enter values:
% colors = {[.8 .25 .25] [.8 .5 .25] [.4 .5 .8] [.25 .25 .9]};
% [h1, s1] = barplot_colored(dat, 'x', [1 2 4 5], 'colors', colors);
%
% Set X Tick Label:
% [h1, s1] = barplot_colored(dat, 'XTicklabels', {'A' 'B' 'C' 'D'});

if iscell(data)
    for k=1:length(data)
        means(k)=mean(data{k});
        stderr(k)=std(data{k})/sqrt(length(data{k}));
    end
else
    means=mean(data);
    stderr=std(data)/sqrt(size(data,1));
end

x = 1:length(means); % can replace x values
myfontsize = 18;
mytitle = '';

for k=1:2:length(varargin)
    if strcmp(varargin{k},'colormap')
        eval(['colormapfun=@' varargin{k+1} ';'])
    elseif strcmp(varargin{k},'title')
        mytitle=varargin{k+1};
    elseif strcmp(lower(varargin{k}),'fontsize')
        myfontsize=varargin{k+1};
    elseif strcmp(lower(varargin{k}),'xticklabel') || strcmp(lower(varargin{k}),'xticklabels')
        XTickLabel=varargin{k+1};
    elseif strcmp(lower(varargin{k}),'ylabel')
        Ylabel=varargin{k+1};
    elseif strcmp(lower(varargin{k}),'xlabel')
        Xlabel=varargin{k+1};
   elseif strcmp(lower(varargin{k}),'title')
        mytitle=varargin{k+1};
    elseif strcmp(varargin{k},'x')
        x=varargin{k+1};
     elseif strcmp(varargin{k},'colors')
        colors=varargin{k+1};
        
    elseif strcmp(varargin{k},'within')
        if iscell(data), error('Within error bars not implemented for cell input data'); end
        
        stderr = barplot_get_within_ste(data);
        stderr = repmat(stderr, 1, length(means));
    end
end

k = size(x, 2);

if ~exist('colormapfun','var')
    colormapfun=@hsv;
end

if ~exist('colors','var')
    colors = colormapfun(k);
end

if ~iscell(colors)
    for i = 1:size(colors, 1)
        tmp{i} = colors(i, :);
    end
    colors = tmp;
end

axh = gca;
set(axh, 'FontSize', myfontsize)
title(mytitle);

for i = 1:k
    hold on
    h(i) = bar(x(i), means(i));
    set(h(i), 'FaceColor', colors{i});
end

% This is old code, Matlab has changed...2014/5 update
% h=bar(x, means);
% s=get(h,'Children');
% 
% colormap(colormapfun(length(means)));
% set(s,'CData',1:length(means));


errorbar(x, means,stderr,'k','LineWidth',2,'LineStyle','none')
set(gca,'Xlim',[0 max(x)+1])
if exist('XTickLabel','var')
    set(gca,'XTickLabel',XTickLabel, 'XTick', x)
else
    set(gca,'XTickLabel',[])
end
if exist('mytitle','var')
    title(mytitle)
end
if exist('Ylabel','var')
    ylabel(Ylabel)
end
if exist('Xlabel','var')
    xlabel(Xlabel)
end

hold off