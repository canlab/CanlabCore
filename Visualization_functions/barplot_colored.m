function [h,s]=barplot_colored(data,varargin)

% haven't had time to document this yet
% this is a good function though
% within-subject error bars now added; use 'within'
%
% [h,s]=barplot_colored(data,varargin)
% 
% h is the axes handle, s is the handle to its children
% 
% Var args: strings, followed by values for each
%     if strcmp(varargin{k},'colormap')
%         eval(['colormapfun=@' varargin{k+1} ';'])
%     elseif strcmp(varargin{k},'title')
%         title=varargin{k+1};
%     elseif strcmp(varargin{k},'XTickLabels')
%         XTickLabels=varargin{k+1};
%     elseif strcmp(varargin{k},'ylabel')
%         ylabel=varargin{k+1};
%     elseif strcmp(varargin{k},'xlabel')
%         xlabel=varargin{k+1};
%     elseif strcmp(varargin{k},'x')
%         x=varargin{k+1};
%      elseif strcmp(varargin{k},'within') -> DO within-subject STE
%     end
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
% colormap([1 0 0; 0 0 1; 1 0 0; 0 0 1])

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

for k=1:2:length(varargin)
    if strcmp(varargin{k},'colormap')
        eval(['colormapfun=@' varargin{k+1} ';'])
    elseif strcmp(varargin{k},'title')
        mytitle=varargin{k+1};
    elseif strcmp(varargin{k},'XTickLabel') || strcmp(varargin{k},'XTickLabels')
        XTickLabel=varargin{k+1};
    elseif strcmp(varargin{k},'Ylabel')
        Ylabel=varargin{k+1};
    elseif strcmp(varargin{k},'Xlabel')
        Xlabel=varargin{k+1};
    elseif strcmp(varargin{k},'x')
        x=varargin{k+1};
        
    elseif strcmp(varargin{k},'within')
        if iscell(data), error('Within error bars not implemented for cell input data'); end
        
        stderr = barplot_get_within_ste(data);
        stderr = repmat(stderr, 1, length(means));
    end
end


if ~exist('colormapfun','var')
    colormapfun=@hsv;
end

h=bar(x, means);
s=get(h,'Children');

colormap(colormapfun(length(means)));
set(s,'CData',1:length(means));

hold on

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