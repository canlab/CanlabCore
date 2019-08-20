function [handle, xdata_all] = bar_wani_2016(y, e, bar_width, varargin)
% Draw a bar plot with error bars with some additional useful features 
% (doesn't work with older than matlab 2015. If you're using matlab 2014 
%  or older, please use bar_wani, not this).  
%
% :Usage:
% ::
%
%    h = bar_wani(y, e, bar_width, varargin)
%
% :Inputs:
%
%   **y:**
%        y values for bars (row: bar grouping, column: bar values 
%        within a bar group) (e.g., if there are m bar groups and 
%        n bars for each group, y should be a m x n matrix)
%
%   **e:**
%        error bars (m x n matrix)
%
%   **bar_width:**
%        value for bar width between 0 and 1. This will determine
%        the intervals between bars. 
%
% :Optional Inputs: Enter keyword followed by variable with values
%
%   **'ylim':**
%        y axis range, (e.g., 'ylim', [-1 1])
%
%   **'ytick':**
%        y tick values (e.g., 'ytick', -.002:.001:.002)
%
%   **'errbar_width':**
%        the horizontal width of error bars (e.g., 'errbar_width', 0)
%
%   **'colors':**
%        bar colors: each row determines the color for each bar in order (n x 3 matrix)
%
%   **'ast':**
%        put asterisks according to p values, which should be
%        given. (e.g., 'ast', p [m x n]) *p<.05, **p<.01, ***p<.001
%
%   **'btwlines':**
%        this option puts lines between bar groups. This is
%        followed by the line style (e.g., 'btwlines', '--');
%
%   **'dosave':**
%        followed by savename (e.g., 'dosave', savename)
%
% :Some advanced options:
%
%   **'scatter':**
%        show individual data points, which should be in cell array
%
%   **'text':**
%        this will put a number for each bar (e.g., 'text', round(y) [m x n])
%
%   **'ast_adj_y_pos':**
%        When the asterisk locations (on y axis) for bars with positive
%        values are off, you can adjust it using this option
%
%   **'ast_adj_y_neg':**
%        When the asterisk locations (on y axis) for bars with negative
%        values are off, you can adjust it using this option
%
%   **'ast_adj_x':**
%        When the asterisk locations (on x axis) are off, you 
%        can adjust it using this option
%
%   **'bar_edgecol':**
%        You can use different colors for bar edges (col [n x 3 matrix])
%
%   **'bar_edgewidth':**
%        You can change linewidths for bar edges 
%
% :Output:
%
%   **h:**
%        graphic handles for a bar plot
%
% :Examples: you can see the output in 
% http://wagerlab.colorado.edu/wiki/doku.php/help/core/figure_gallery
% :Examples:
% ::
%
%    % data
%    y = [-0.6518   -0.6934   -0.5417   -0.6496   -0.5946   -0.3839
%        1.1511    0.9090    1.1681    1.2892    0.9346    1.1383];
%    e = [0.3226    0.2936    0.3080    0.3203    0.3368    0.3167
%        0.4026    0.4088    0.4012    0.5586    0.3734    0.4257];
%    p = [0.0433    0.0182    0.0785    0.0426    0.0775    0.2255
%        0.0042    0.0262    0.0036    0.0210    0.0123    0.0075];
% 
%    col =  [0    0.1157    0.2686
%            0.1157    0.2765    0.4725
%            0.4843    0.1157    0.1078
%            0.3667    0.4765    0.1353
%            0.2765    0.1902    0.3824
%            0.0922    0.4216    0.5118
%            0.7941    0.3235   0];
%
%    % draw
%    bar_wani(y, e, .8, 'colors', col, 'errbar_width', 0, 'ast', p, 'ylim', [-2.5 2.5], 'ytick', -2:2, 'ast_adj_x', 0, 'ast_adj_y_neg', .15);
%    set(gca, 'ytickLabel', num2str(get(gca, 'ytick')'));
%    set(gcf, 'position', [1   531   399   169]);
% 
%    savename = 'example_barwani.pdf';
% 
%    try
%        pagesetup(gcf);
%        saveas(gcf, savename);
%    catch
%        pagesetup(gcf);
%        saveas(gcf, savename);   
%    end
%
% ..
%    Copyright (C) 2014  Wani Woo
%
%    Programmers' notes:
%    'yline'
% ..

dosave = 0;
doman_ylim = 0;
docolor = 0;
doast = 0;
doerrbar_width = 0;
doytick = 0;
dobtwlines = 0;
dotext = 0;
doscatter = 0;
ast_adj_y_pos = .2;
ast_adj_y_neg = .5;
ast_adj_x = 0; 
text_adj_y = 0;
bar_edgewidth = 1.8;
doyline = 0;
use_samefig = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'dosave', 'save'}
                dosave = 1;
                savename = varargin{i+1}; 
            case {'ylim'}
                doman_ylim = 1;
                ylim = varargin{i+1};
            case {'colors', 'color'}
                docolor = 1; colors = varargin{i+1};
            case {'ast'}
                doast = 1; p = varargin{i+1};
            case {'errbar_width'}
                doerrbar_width = 1; errwidth = varargin{i+1};
            case {'ytick'}
                doytick = 1; ytick = varargin{i+1};
            case {'btwlines'}
                dobtwlines = 1; btwlnstyle = varargin{i+1};
            case {'yline'}
                doyline = 1; yln_y = varargin{i+1};
            case {'text'}
                dotext = 1; text_bar = varargin{i+1};
            case {'scatter'}
                doscatter = 1; sc_data = fliplr(varargin{i+1});
            case {'ast_adj_y_pos'}
                ast_adj_y_pos = ast_adj_y_pos + varargin{i+1};
            case {'ast_adj_y_neg'}
                ast_adj_y_neg = ast_adj_y_neg + varargin{i+1};
            case {'ast_adj_x'}
                ast_adj_x = ast_adj_x + varargin{i+1};
            case {'text_adj_y'}
                text_adj_y = varargin{i+1};
            case {'bar_edgecol'}
                bar_edgecol = flipud(varargin{i+1});
            case {'bar_edgewidth'}
                bar_edgewidth = varargin{i+1};
            case {'use_samefig'}
                use_samefig = true;
        end
    end
end

if ~use_samefig
    handle.main = create_figure('Bar plot');
    set(gcf, 'Position', [1   450   300.*size(y,1)   250]);
end

barnum = size(y,2);
grnum = size(y,1);

if ~doman_ylim
    step = (max(max(y+e))-min(min(y+e)))*.2;
    ymax = max(max(y+e)) + step;
    ymin = min(min(y-e)) - step;
else
    ymin = ylim(1);
    ymax = ylim(2);
end

if doyline
    line([0.55 .45+size(y,1)], [yln_y yln_y], 'linewidth', 1, 'linestyle', '--', 'color', [.4 .4 .4]);  
end

handle.bar = barweb(y, e, bar_width, [], [], [], [], [], [], [], [], []);
handle.bar = barweb(y, e, bar_width, [], [], [], [], [], [], [], [], []);
set(gcf, 'Color', 'w')
set(gca, 'ylim', [ymin ymax], 'XLim', [0.55 .45+size(y,1)], 'fontsize', 20, 'linewidth', 1.8); % **ADJUST**: adjust basic setting for axis

if doytick
    set(gca, 'ytick', ytick);
end

h = get(gca, 'children');

if docolor
    for i = (barnum+1):2*barnum
        try
            colors_flip = fliplr(colors);
            set(h(i), 'FaceColor', colors_flip{i-barnum});
        catch
            colors_flip = flipud(colors);
            set(h(i), 'FaceColor', colors_flip(i-barnum,:));
        end
    end
end

if doast
    ast = p_asterisk(p);
    for i = 1:barnum
        for j = 1:grnum
            if y(j,i) > 0
                ast_loc(j,i) = y(j,i) + e(j,i) + ast_adj_y_pos * max(max(e)); 
            else
                ast_loc(j,i) = y(j,i) - e(j,i) - ast_adj_y_neg * max(max(e)); 
            end
        end
    end
end

for i = 1:barnum
    if ~doerrbar_width
        errwidth = 8;
    end
    
    handle.bar.errors(i).CapSize = errwidth;
    xdata_all(:,i) = handle.bar.errors(i).XData';
    
    for j = 1:grnum %numel(handle.bar.errors(i).XData)
        if doast
            text(double(handle.bar.errors(i).XData(j)+ast_adj_x), double(ast_loc(j,i)), ast{j,i}, 'fontsize', 18, 'HorizontalAlignment','center'); % **ADJUST** font size
            if dotext
                text_y = ymax*1.15 + text_adj_y;
                text(handle.bar.errors(i).XData(j), text_y, num2str(text_bar(j,i)), 'fontsize', 20, 'HorizontalAlignment','center'); % **ADJUST** font size
            end
        end
        
        if doscatter
            hold on; scatter(double((handle.bar.errors(i).XData(j))+.02).*ones(size(sc_data{j,i})), sc_data{j,i}, 50, [.3 .3 .3], 'filled');
            if ymin > min(sc_data{j,i}), ymin = min(sc_data{j,i}); set(gca, 'ylim', [ymin ymax]); end
            if ymax < max(sc_data{j,i}), ymax = max(sc_data{j,i}); set(gca, 'ylim', [ymin ymax]); end
        end
    end
    
end

j = 1; 
if ~exist('bar_edgecol', 'var'), bar_edgecol = repmat([0 0 0], barnum, 1); end
for i = barnum+1:2*barnum, set(h(i), 'lineWidth', bar_edgewidth, 'edgecolor', bar_edgecol(j,:)); j = j+1; end % **ADJUST** bar linewidth and edge color

% btwlines: drawing lines between bar groups
if dobtwlines
    btwx = 1.5:1:grnum;
    btwy = get(gca, 'ylim');
    btwx = repmat(btwx', 1, 2);
    btwy = repmat(btwy, size(btwx,1),1);
    btwlnwidth = 1;
    btwcol = [.3 .3 .3];
    for i = 1:size(btwx,1)
        line(btwx(i,:), btwy(i,:), 'linestyle', btwlnstyle, 'linewidth', btwlnwidth, 'color', btwcol)
    end
end

set(gca, 'tickdir', 'out', 'ticklength', [.01 .01]);

xtickdata = sort(xdata_all(:));

if size(y,1)==1
    set(gca, 'xtick', xtickdata(~isnan(y)));
else
    set(gca, 'xtick', xtickdata(~any(isnan(y(:)),2)));
end

if dosave
    try
        pagesetup(handle);
        saveas(handle, savename);
    catch
        pagesetup(handle);
        saveas(handle, savename);
    end 
end

xdata_all = flipud(xdata_all);

end

function ast = p_asterisk(p)
% function ast = p_asterisk(p)
% change p values into asterisk
% 

p = double(p < .1) + double(p < .05) + double(p < .01) + double(p < .001);

ast = cell(size(p));

for i = 1:size(p,1)
    for j = 1:size(p,2)
        if p(i,j) > 1
            ast{i,j} = repmat('*', 1, p(i,j)-1);
        elseif p(i,j) == 1
            ast{i,j} = '';
        else
            ast{i,j} = '';
        end
    end
end

end


function handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)

% This function is from http://www.mathworks.com/matlabcentral/fileexchange/10803-barweb--bargraph-with-error-bars-

% Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
%
% Ex: handles = barweb(my_barvalues, my_errors, [], [], [], [], [], bone, [], bw_legend, 1, 'axis')
%
% barweb is the m-by-n matrix of barvalues to be plotted.
% barweb calls the MATLAB bar function and plots m groups of n bars using the width and bw_colormap parameters.
% If you want all the bars to be the same color, then set bw_colormap equal to the RBG matrix value ie. (bw_colormap = [1 0 0] for all red bars)
% barweb then calls the MATLAB errorbar function to draw barvalues with error bars of length error.
% groupnames is an m-length cellstr vector of groupnames (i.e. groupnames = {'group 1'; 'group 2'}).  For no groupnames, enter [] or {}
% The errors matrix is of the same form of the barvalues matrix, namely m group of n errors.
% Gridstatus is either 'x','xy', 'y', or 'none' for no grid.
% No legend will be shown if the legend paramter is not provided
% 'error_sides = 2' plots +/- std while 'error_sides = 1' plots just + std
% legend_type = 'axis' produces the legend along the x-axis while legend_type = 'plot' produces the standard legend.  See figure for more details
%
% The following default values are used if parameters are left out or skipped by using [].
% width = 1 (0 < width < 1; widths greater than 1 will produce overlapping bars)
% groupnames = '1', '2', ... number_of_groups
% bw_title, bw_xlabel, bw_ylabel = []
% bw_color_map = jet
% gridstatus = 'none'
% bw_legend = []
% error_sides = 2;
% legend_type = 'plot';
%
% A list of handles are returned so that the user can change the properties of the plot
% handles.ax: handle to current axis
% handles.bars: handle to bar plot
% handles.errors: a vector of handles to the error plots, with each handle corresponding to a column in the error matrix
% handles.legend: handle to legend
%
% See the MATLAB functions bar and errorbar for more information
%
% Author: Bolu Ajiboye
% Created: October 18, 2005 (ver 1.0)
% Updated: Dec 07, 2006 (ver 2.1)
% Updated: July 21, 2008 (ver 2.3)

% Get function arguments
if nargin < 2
	error('Must have at least the first two arguments:  barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, barwebtype)');
elseif nargin == 2
	width = 1;
	groupnames = 1:size(barvalues,1);
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 3
	groupnames = 1:size(barvalues,1);
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 4
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 5
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 6
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 7
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 8
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 9
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 10
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == 11
	legend_type = 'plot';
end

change_axis = 0;
ymax = 0;

if size(barvalues,1) ~= size(errors,1) || size(barvalues,2) ~= size(errors,2)
	error('barvalues and errors matrix must be of same dimension');
else
	if size(barvalues,2) == 1
		barvalues = barvalues';
		errors = errors';
	end
	if size(barvalues,1) == 1
		barvalues = [barvalues; zeros(1,length(barvalues))];
		errors = [errors; zeros(1,size(barvalues,2))];
		change_axis = 1;
	end
	numgroups = size(barvalues, 1); % number of groups
	numbars = size(barvalues, 2); % number of bars in a group
	if isempty(width)
		width = 1;
	end
	
	% Plot bars
	handles.bars = bar(barvalues, width,'edgecolor','k', 'linewidth', 2);
    h_temp = handles.bars(any(isnan(barvalues)));
    set(h_temp, 'linestyle', 'none')
    
	hold on
	if ~isempty(bw_colormap)
		colormap(bw_colormap);
	else
		colormap(jet);
	end
	if ~isempty(bw_legend) && ~strcmp(legend_type, 'axis')
		handles.legend = legend(bw_legend, 'location', 'best', 'fontsize',12);
		legend boxoff;
	else
		handles.legend = [];
	end
	
	% Plot erros
    % 	for i = 1:numbars
    % 		x =get(get(handles.bars(i),'children'), 'xdata');
    % 		x = mean(x([1 3],:));
    % 		handles.errors(i) = errorbar(x, barvalues(:,i), errors(:,i), 'k', 'linestyle', 'none', 'linewidth', 2);
    % 		ymax = max([ymax; barvalues(:,i)+errors(:,i)]);
    %         ymin = min([ymax; barvalues(:,i)+errors(:,i)]); % wani added
    % 	end
    for i = 1:numbars
        if ~verLessThan('matlab', '8.4') % HG2
            x =  handles.bars(i).XData + handles.bars(i).XOffset;
        else
            x =get(get(handles.bars(i),'children'), 'xdata');
            x = mean(x([1 3],:));
        end
        handles.errors(i) = errorbar(x, barvalues(:,i), errors(:,i), 'k', 'linestyle', 'none', 'linewidth', 2);
        ymax = max([ymax; barvalues(:,i)+errors(:,i)]);
        ymin = min([ymax; barvalues(:,i)+errors(:,i)]); % wani added
    end
	
    if error_sides == 1
        set(gca,'children', flipud(get(gca,'children')));
    end
    
    ylim([ymin*1.1 ymax*1.1]); % wani changed
    
	xlim([0.5 numgroups-change_axis+0.5]);
	
	if strcmp(legend_type, 'axis')
		for i = 1:numbars
			xdata = get(handles.errors(i),'xdata');
			for j = 1:length(xdata)
				text(xdata(j),  -0.03*ymax*1.1, bw_legend(i), 'Rotation', 60, 'fontsize', 12, 'HorizontalAlignment', 'right');
			end
		end
		set(gca,'xaxislocation','top');
	end
	
	if ~isempty(bw_title)
		title(bw_title, 'fontsize',14);
	end
	if ~isempty(bw_xlabel)
		xlabel(bw_xlabel, 'fontsize',14);
	end
	if ~isempty(bw_ylabel)
		ylabel(bw_ylabel, 'fontsize',14);
	end
	
	set(gca, 'xticklabel', groupnames, 'box', 'off', 'ticklength', [0 0], 'fontsize', 12, 'xtick',1:numgroups, 'linewidth', 2,'xgrid','off','ygrid','off');
	if ~isempty(gridstatus) && any(gridstatus == 'x')
		set(gca,'xgrid','on');
	end
	if ~isempty(gridstatus) && any(gridstatus ==  'y')
		set(gca,'ygrid','on');
	end
	
	handles.ax = gca;
	
	hold off
end

end
