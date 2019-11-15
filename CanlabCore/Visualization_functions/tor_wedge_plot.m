function [handles, key_points] = tor_wedge_plot(radius_values, text_labels, varargin)
% Make polar wedge plot
% - Single column of k radius_values makes plot with k wedges
% - Matrix of n obs x k variables creates plot with mean +/- se for k wedges
% - Negative values: Default is to plot stripes for neg values
% - Option to make them a different color
% - Option to not show negative relationships by setting negative radius values to zero
% - image_similarity_plot.m uses tor_wedge_plot to plot pos and neg associations in different colors
%   (see example below)
%
% :Usage:
% ::
%
%   [handles, key_points] = tor_wedge_plot(radius_values, colors, text_labels, [optional arguments])
%
% :Inputs:
%
%   **radius_values:**
%        vector of radius values
%
%   **names:**
%        is cell array, one cell per wedge, with names for each condition
%
% :Optional Inputs:
%
%   **'nofigure':**
%        suppress figure
%
%   **'labelstyle':**
%        followed by string: 'equal' (default) or 'close'
%
%   **'nocircle':**
%        suppress outer reference circle
%
%   **'nospokes':**
%        suppress reference spokes
%
%   **'ignorenegative':**
%        true or false [default]. Ignore negative values by setting them to
%        zero
%
%   **'bicolor':**
%        two-color plot; with positive-valued entries in one color and
%        negative-valued ones in another. Enter 2 colors in colors cell (see below).
%        Will plot positive values in color{1} and negative values in color{2}
%
%   **'colors':**
%        Followed by cell array, one cell per wedge, with [r g b] color triplets
%        OR
%        [{pos_color_triplet} {neg_color_triplet}] when using 'bicolor' option
%
%
%   **'colorband_colors':**
%        Followed by cell array with one nested cell per wedge, with [r g b] color triplets
%
%
%   **'linewidth':**
%       Followed by line width
%
%   **'outer_circle_radius':**
%       Followed by radius for outer guide circle
%
%   **'nofill':**
%       Turn off color fill in wedges
%
%   **'errorbars'**
%       True or false; Default for matrix input
%       Plot wedges for mean +/- standard error
%       Treats columns as variables for wedges, rows as replicates
%       Similar function as 'average' in image_similarity_plot.m
%
% :Output:
%
%   **handles:**
%        Handles to line and fill objects; see draw_pie_wedge.m
%
%   **key_points:**
%        Handles to fill objects; see draw_pie_wedge.m
%
% :Examples:
% ::
%
% % Create a made-up plot of similarity with brain networks:
% --------------------------------------------------------------
% [obj, netnames, imgnames] = load_image_set('bucknerlab');
% text_labels = netnames;
% radius_values = (1:7) ./ 7;
% [handles, key_points] = tor_wedge_plot(radius_values, text_labels);
%
% % Change colors and options:
% [handles, key_points] = tor_wedge_plot(radius_values, text_labels, 'colors', scn_standard_colors(7), 'labelstyle', 'close', 'nospokes', 'nocircle');
%
% Example of how to show positive and negative associations in different
% colors on the same plot (see image_similarity_plot.m)
% --------------------------------------------------------------
% mycolors = repmat({[1 .7 0]}, 1, length(networknames));
% wh = r < 0;
%
% % plot negative values in the complementary color
% mycolors(wh) = repmat({[1 1 1] - groupColors{1}}, 1, sum(wh));
% r_to_plot = abs(r);
%
% outercircleradius = .3;
%
% hh = tor_wedge_plot(r_to_plot, networknames, 'outer_circle_radius', outercircleradius, 'colors', mycolors, 'nofigure');



% Default values for options
% -------------------------------------------------------------------------
input_args = varargin;                  % save for later to pass in recursively
ignorenegative = false;                 % ignore negative-valued elements by zeroing them out
bicolor = false;                        % two-color plot, with positive-valued entries in one color and negative-valued ones in another (enter 2 colors in colors cell)
dofigure = true;
docircle = true;
dospokes = true;
labelstyle = 'equal';  % or 'close'
linewidth = 1;
outer_circle_radius = 1;
nofill = false;
docolorband=false;

radius_values = double(radius_values);

if min(size(radius_values)) > 1
    n_categories = size(radius_values, 2);
else
    n_categories = length(radius_values);
end

dostripes = false(1, n_categories);     % striping/shading, default for negative values

try
    colors = seaborn_colors(n_categories + 1); % add one to make colors more diff
catch %in case we have many maps
    colors = hsv(n_categories + 1); % add one to make colors more diff
end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'ignorenegative', ignorenegative = true;
                
            case 'bicolor', bicolor = true;
                
            case 'nofigure', dofigure = false;
                
            case 'nocircle', docircle = false;
                
            case 'nospokes', dospokes = false;
                
            case 'nofill', nofill = true;
                
            case 'colorband', docolorband = true;
                
            case 'labelstyle', labelstyle = varargin{i+1}; varargin{i+1} = [];
                
            case 'colors', colors = varargin{i+1}; varargin{i+1} = [];
                
            case 'linewidth', linewidth = varargin{i+1}; varargin{i+1} = [];
                
            case 'outer_circle_radius', outer_circle_radius = varargin{i+1}; varargin{i+1} = [];
                
            case 'colorband_colors', colorband_colors = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if ~docolorband
    colorband_colors = scn_standard_colors(n_categories);
end
    
% Error checks

if any(size(outer_circle_radius)) > 1, error('''outer_circle_radius'' input should be followed by a scalar value.'); end


% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% Recursive function callback
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------

% Handle matrix input (replicates = rows) by calling this function
% recursively. Makes plot of mean + error zones
% --------------------------------------------------------------------------

if min(size(radius_values)) > 1  % is matrix
    
    clear handles key_points
    
    % transpose if only the number of rows == num labels.  Columns should match with labels.
    mysz = size(radius_values);
    
    if mysz(1) == length(text_labels) && mysz(2) ~= length(text_labels)
        % could be transposed; try transposing
        radius_values = radius_values';
    end
    
    m = nanmean(radius_values)';
    se = ste(radius_values)';
    
    % manually build colors to handle both pos and neg values
    isneg = m < 0;
    mycolors = repmat(colors(1), 1, n_categories);
    mycolors(isneg) = colors(2);
    
    % get rid of bicolor input because we've manually specified colors
    wh = strcmp(input_args, 'bicolor'); input_args(wh) = [];
    
    % outer wedge + se
    [handles.outer, key_points.outer] = tor_wedge_plot(abs(m)+se, text_labels, input_args{:}, 'colors', mycolors);
    
    % inner bound, - se
    [handles.inner, key_points.inner] = tor_wedge_plot(abs(m)-se, text_labels, 'nocircle', 'nofigure', input_args{:}, 'colors', {[1 1 1] [1 1 1]}, 'ignorenegative'); % override varargin with later inputs; need 2 colors for wedge
    set([handles.inner(:).fill_han], 'FaceAlpha', .8);
    
    % Mean line, no fill
    [handles.meanline, key_points.meanline] = tor_wedge_plot(abs(m), text_labels, 'nofigure', input_args{:}, 'linewidth', 3, 'nocircle', 'nofill', 'colors', mycolors);
    
    % delete redundant text
    delete([handles(:).inner.texth])
    [handles(:).inner.texth] = deal([]);
    delete([handles(:).meanline.texth])
    [handles(:).meanline.texth] = deal([]);
    
    mymeanlines = cat(1, handles.meanline.line_han);
    delete(mymeanlines(:, 2)')
    delete(mymeanlines(:, 3)')
    for i = 1:length(handles.meanline), handles.meanline(i).line_han(2:3) = []; end
    
    return
end

% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% End recursive function callback - Main function below
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------

% Handle color, stripes, and neg entry options
% --------------------------------------------------------------------------

if bicolor
    % plot positive values in color{1} and negative values in color{2}
    
    mycolors = repmat(colors(1), 1, n_categories);
    wh = radius_values < 0;
    
    % plot negative values in the complementary color
    mycolors(wh) = colors(2);  % repmat({[1 1 1] - groupColors{1}}, 1, sum(wh));
    
    colors = mycolors;
    
else  % stripes for negative relationships
    
    if length(colors) < n_categories
        colors = repmat(colors(1), 1, n_categories);
    end
    
    % Put stripes on negative-valued entries
    dostripes(radius_values < 0) = true;
end

if ignorenegative
    
    % Will not show negative relationships; zero them out
    radius_values(radius_values < 0) = 0;
    
else
    % Make all values positive, because we have either color- or
    % stripe-coded them
    
    radius_values = abs(radius_values);
end

% Calculations

st = linspace(0, 2*pi, n_categories+1); % starting vals in radians

en = st(2:end);
st = st(1:end-1);

clear handles key_points

if dofigure
    create_figure('wedge'); hold on;
end

% Reference circle
% -------------------------------------------------------------------------
if docircle
    
    circle_handles = draw_pie_wedge(0, 2*pi, outer_circle_radius, 'linecolor', [.5 .5 .5], 'fillcolor', 'none');
    
    delete(circle_handles.line_han(2:3))
    circle_handles.line_han = circle_handles.line_han(1);
    set(circle_handles.line_han, 'LineWidth', .75);
    
end

% Spokes
% -------------------------------------------------------------------------
if dospokes
    
    for i = 1:n_categories
        
        spoke_handles = draw_pie_wedge(st(i), en(i), outer_circle_radius, 'linecolor', [.5 .5 .5], 'fillcolor', 'none');
        
        delete(spoke_handles.line_han([1 3]))
        spoke_handles.line_han = spoke_handles.line_han(2);
        set(spoke_handles.line_han, 'LineWidth', .75, 'LineStyle', ':');
    end
    
end

% Wedges
% -------------------------------------------------------------------------

for i = 1:n_categories
    
    if nofill
        
        [handles(i), key_points(i)] = draw_pie_wedge(st(i), en(i), radius_values(i), ...
            'linecolor', colors{i} ./ 2, 'fillcolor', 'none', 'linewidth', linewidth, ...
            'stripes', dostripes(i), 'stripedensity', 30);
        
    else
        
        [handles(i), key_points(i)] = draw_pie_wedge(st(i), en(i), radius_values(i), ...
            'linecolor', colors{i} ./ 2, 'fillcolor', colors{i}, 'linewidth', linewidth, ...
            'stripes', dostripes(i), 'stripedensity', 30);
        
    end
end

axis equal

% Add text labels
% -------------------------------------------------------------------------
% Get rotation from tangent angle
myrotation = [key_points(:).tangentangle];
wh = myrotation > 90 & myrotation < 270;
myrotation(wh) = myrotation(wh) - 180;

handles(i).texth = [];

%mytextsize = 16 * 6/n_categories;
mytextsize = 30 ./ (n_categories.^.3);

switch labelstyle
    
    case 'close'
        
        for i = 1:n_categories
            handles(i).texth = text(key_points(i).xmid_outside, key_points(i).ymid_outside, text_labels{i}, ...
                'Rotation', myrotation(i), 'FontSize', mytextsize,...
                'HorizontalAlignment','center','parent',gca, 'Layer','front');
        end
        
    case 'equal'
        
        
        for i = 1:n_categories
            tmid = key_points(i).tmid_radians;
            
            xmid_outside = 0 + (outer_circle_radius + .04*outer_circle_radius) * cos(tmid);
            ymid_outside = 0 + (outer_circle_radius + .04*outer_circle_radius) * sin(tmid);
            
            handles(i).texth = text(xmid_outside, ymid_outside, text_labels{i}, ...
                'Rotation', myrotation(i), 'FontSize', mytextsize,...
                'HorizontalAlignment','center','parent',gca, 'Layer','front');
            
        end
        
    case 'curvy'
        breakpoints=0:2*pi/n_categories:2*pi;

        for i = 1:n_categories
            
            text_location_handles{i} = draw_pie_wedge(breakpoints(i), breakpoints(i+1), outer_circle_radius+ .28*outer_circle_radius, 'linecolor', 'none', 'fillcolor', 'none');
            delete(text_location_handles{i}.line_han(2:3))
            
            xy = fliplr([text_location_handles{i}.line_han(1).XData;
                text_location_handles{i}.line_han(1).YData]);
            if size(xy,1)>2; xy=xy'; end
            
            xy = double(xy); % tor: added to avoid errors in text()
            
            m = length(text_labels{i});
            xy = xy(:,floor(max(1,(size(xy,2)/2-4*m))):floor(min((size(xy,2)/2+4*m),size(xy,2)))); %squeeze together a bit
            
            n = size(xy,2);
            
            XY = spline(1:n,xy,linspace(1,n,m+1));
            dXY = XY(:,2:end)-XY(:,1:end-1);
            theta = (arrayfun(@(y,x) atan2(y,x),dXY(2,:),dXY(1,:)))/2/pi*360;
            
            XY = (XY(:,1:end-1)+XY(:,2:end))/2;
            hold on;
            for ii=1:m
                handles(i).texth(ii)=text(XY(1,ii),XY(2,ii),text_labels{i}(ii),'rotation',theta(ii),...
                    'horizontalalignment','center','verticalalignment','bottom', 'FontSize', mytextsize);
            end
        
        end
    case 'radial'
       breakpoints=0:2*pi/n_categories:2*pi;

        for i = 1:n_categories
            
            text_location_handles{i} = draw_pie_wedge(breakpoints(i), breakpoints(i+1), outer_circle_radius+ .28*outer_circle_radius, 'linecolor', 'none', 'fillcolor', 'none');
            delete(text_location_handles{i}.line_han(2:3))
            
            xy = fliplr([text_location_handles{i}.line_han(1).XData;
                text_location_handles{i}.line_han(1).YData]);
            if size(xy,1)>2; xy=xy'; end
            
            xy = double(xy); % tor: added to avoid errors in text()
            
            m = length(text_labels{i});
            xy = xy(:,floor(max(1,(size(xy,2)/2-4*m))):floor(min((size(xy,2)/2+4*m),size(xy,2)))); %squeeze together a bit
            
            n = size(xy,2);
            
            XY = spline(1:n,xy,linspace(1,n,m+1));
            dXY = XY(:,2:end)-XY(:,1:end-1);
            theta = (arrayfun(@(y,x) atan2(y,x),dXY(2,:),dXY(1,:)))/2/pi*360 + 90;
            
            XY = (XY(:,1:end-1)+XY(:,2:end))/2;
            hold on;
            
            nn = length(XY);
            midpoint = round(nn/2);
            if theta(midpoint) > 90
                theta = theta(midpoint) - 180;
                alignment = 'right';
            else
                theta = theta(midpoint);
                alignment = 'left';
            end
            handels(i).texth(1) = text(XY(1,midpoint),XY(2,midpoint),text_labels{i},'rotation',theta,...
                'horizontalalignment',alignment,'verticalalignment','bottom', 'FontSize', mytextsize);
        
        end         
        
    otherwise
        error('Unknown labelstyle');
        
end % labelstyle

axis off

% Add colorband
% -------------------------------------------------------------------------

if docolorband
    
    % for each wedge
    breakpoints=0:2*pi/n_categories:2*pi;
    
    if length(colorband_colors) ~= n_categories
    
        disp('Warning: tor_wedge_plot: Number of colorband_colors does not equal number of categories to plot.');
        fprintf('Length colorband_colors = %3.0f, n_categories = %3.0f\nSkipping color band. Check inputs or code.\n', length(colorband_colors), n_categories);
        
    else 
        
        for i=1:n_categories
            
            colorband_handles{i} = draw_pie_wedge(breakpoints(i), breakpoints(i+1), outer_circle_radius+ .2 * outer_circle_radius, 'linecolor', colorband_colors{i}, 'fillcolor', 'none');
            
            delete(colorband_handles{i}.line_han(2:3))
            colorband_handles{i}.line_han = colorband_handles{i}.line_han(1);
            set(colorband_handles{i}.line_han, 'LineWidth', 10);
            
        end
        
    end
    
end

end % main function