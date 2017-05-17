function [f1, axh] = create_figure(tagname, varargin)
% :Usage:
% ::
%
%    f1 = create_figure(['tagname'], [subplotrows], [subplotcols], [do not clear flag])
%
% checks for old figure with tag of tagname,
% clears it if it exists, or creates new one if it doesn't

% Edit: 1/2017 by Tor Wager - set figure aspect ratio in proportion to data
% plots, when creating new figure

axh = [];

if nargin < 1 || isempty(tagname)
    tagname = 'nmdsfig';
end

doclear = 1;    % clear if new or if old and existing
createnew = 0;

if length(varargin) > 2 && varargin{3}
    % use same figure; do not clear
    doclear = 0;
end

old = findobj('Tag', tagname);
old = old( strcmp( get(old, 'Type'), 'figure' ) );

if ~isempty(old)
    
    if length(old) > 1
        % multiple figures with same tag!
        close(old(2:end))
        old = old(1);
    end
    
    if doclear, clf(old); end
    
    f1 = old;
    
else
    % Or create new
    createnew = true;
    
    %         scnsize = get(0,'ScreenSize');
    %
    %         xdim = min(scnsize(3)./2, 700);
    %         ydim = min(scnsize(4)./2, 700);
    %
    %         f1 = figure('position',round([50 50 xdim ydim]),'color','white');
    
    f1 = figure; %('color','white');
    
    set(f1, 'Tag', tagname, 'Name', tagname, 'color', 'white');
    hold on
    
end

% activate this figure
figure(f1);


if doclear % true for new figs or cleared ones
    
    % Create subplots, if requested; set axis font sizes
    
    if length(varargin) > 0
        i = max(1, varargin{1});
        j = max(1, varargin{2});
    else
        i = 1;
        j = 1;
    end
    
    np = max(1, i * j);
    
    % quit and return if no subplots
    if np == 1
        axh = gca;
        cla
        set(gca,'FontSize', 14)
        hold on
        return
    end
    
    for k = 1:np
        
        axh(k) = subplot(i,j,k);
        cla;
        set(gca,'FontSize', 14)
        hold on
    end
    
    axes(axh(1));
    
end % if doclear

if createnew
    
    set_figure_position(i, j);
    
end

end % function


function  set_figure_position(i, j)
figaspect = j / i;  % aspect ratio, width / height
pos = get(gcf, 'Position');
w = pos(3);
h = pos(4);
w = w .* figaspect;


screensz = get(0, 'screensize');
maxw = screensz(3) .* .8;
maxh = screensz(4) .* .8;

reducefactor = max(w ./ maxw, h ./ maxh);  % in case we go beyond screen boundaries

pos(1:2) = [30 screensz(4) - 30]; % set to top corner

pos(3) = w ./ reducefactor;
pos(4) = h ./ reducefactor;

set(gcf, 'Position', pos);

end
