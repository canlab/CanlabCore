function [h]=scatter_glm(x,y,stats,varargin)
% scatterplots x and y and draws the logistic function specified by the
% parameters in the stats structure (output from glmfit).
%
% :Usage:
% ::
%
%    [h]=scatter_glm(x,y,stats,varargin)
%
% :Output:
%
%   **h:**
%        is the handle to the figure axis.

Title = [];
Xlabel = [];
Ylabel = [];
fontsize = 18;
xfontsize = 20;
yfontsize = 20;
tfontsize = 24;
param=1;



for k = 1:length(varargin)
    if ischar(varargin{k})
        switch(varargin{k})
            case 'Title'
                Title=varargin{k+1};
            case 'XLabel'
                Xlabel=varargin{k+1};
            case 'Ylabel'
                Ylabel=varargin{k+1};
            case 'FontSize'
                fontsize=varargin{k+1};
            case 'xFontSize'
                xfontsize=varargin{k+1};
            case 'yFontSize'
                yfontsize=varargin{k+1};
            case 'tFontSize'
                beep,disp('Title Font size control has not yet been implemented')
                tfontsize=varargin{k+1};
            case 'param'
                param=varargin{k+1};
        end
    end
end


figure;
h=scatter(x,y,'filled','MarkerEdgeColor','k','MarkerFaceColor','k');
h=get(h,'Parent');
xlabel(Xlabel,'FontSize',xfontsize)
ylabel(Ylabel,'FontSize',yfontsize)
hold on
lim=get(h,'XLim');
x=lim(1):(lim(2)-lim(1))/1000:lim(2);
y=exp(stats.beta(1)+stats.beta(param+1)*x)./(1+exp(stats.beta(1)+stats.beta(param+1)*x));
plot(x,y,'k')
hold off
