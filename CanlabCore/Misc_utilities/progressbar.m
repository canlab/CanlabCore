function [f,ax] = progressbar(meth,val)
% Create a progress bar window with tag 'progressbar'
% that can be updated
%
% :Usage:
% ::
%
%     progressbar(meth,val)
%
% :Inputs:
%
%   **meth:**
%        can be 'init' or 'update'
%
% x-axis limits are 0 - 100, so val should be % complete for best results
%
% :Outputs:
%
%   the f and ax handles are not really needed
%   as 'update' finds the axis with 'progressbar' tag.
%
% :Example:
% ::
%
%    progressbar('update',100*i./nvars);
%
% ..
%    tor wager
% ..

    f = [];
    ax = [];

    switch meth
        case 'init'
            f = figure('Color','w');
            tmp = get(gcf,'Position') .* [1 1 .5 .1];
            set(gcf,'Position',tmp)
            set(gcf,'MenuBar','none','NumberTitle','off')
            figure(f), set(gca,'Xlim',[0 100])
            subplot(1,1,1);
            ax = findobj('Type','Axes');
            set(ax,'XLim',[0 100],'YLim',[0 1]);
            set(ax,'Tag','progressbar')
            drawnow
            
        case 'update'
            ax = findobj('Tag','progressbar');
            if isempty(ax), progressbar('init'); end
            delete(findobj(ax,'Type','patch'));  % delete current bars
            barh(val);
            set(ax,'XLim',[0 100]);
            drawnow
            
        otherwise
            error('progressbar: Unknown method');
    end

    return
