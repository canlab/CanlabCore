function [handles, key_points] = tor_wedge_plot(radius_values, text_labels, varargin)
% Make polar wedge plot
% - Will not show negative relationships; zeroes them out
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
%   **'colors':**
%        Followed by cell array, one cell per wedge, with [r g b] color triplets
%
%   **'linewidth':**
%       Followed by line width
%
%   **'outer_circle_radius':**
%       Followed by radius for outer guide circle
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



% Defaults and inputs
dofigure = true;
docircle = true;
dospokes = true;
labelstyle = 'equal';  % or 'close'
n_categories = length(radius_values);
linewidth = 1;
outer_circle_radius = 1;

% Will not show negative relationships; zero them out
radius_values(radius_values < 0) = 0;

colors = seaborn_colors(n_categories + 1); % add one to make colors more diff

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'nofigure', dofigure = false;
                
            case 'nocircle', docircle = false;
                
            case 'nospokes', dospokes = false;
                
            case 'labelstyle', labelstyle = varargin{i+1}; varargin{i+1} = [];
                
            case 'colors', colors = varargin{i+1}; varargin{i+1} = [];
                
            case 'linewidth', linewidth = varargin{i+1}; varargin{i+1} = [];
                
            case 'outer_circle_radius', outer_circle_radius = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
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
    
    [handles(i) key_points(i)] = draw_pie_wedge(st(i), en(i), radius_values(i), 'linecolor', colors{i} ./ 2, 'fillcolor', colors{i}, 'linewidth', linewidth);
    
end

axis equal

% Add text labels
% -------------------------------------------------------------------------
% Get rotation from tangent angle
myrotation = [key_points(:).tangentangle];
wh = myrotation > 90 & myrotation < 270;
myrotation(wh) = myrotation(wh) - 180;

handles(i).texth = [];

switch labelstyle
    
    case 'close'
        
        for i = 1:n_categories
            handles(i).texth = text(key_points(i).xmid_outside, key_points(i).ymid_outside, text_labels{i}, ...
                'Rotation', myrotation(i), 'FontSize', 16 * 6/n_categories,...
                'HorizontalAlignment','center','parent',gca, 'Layer','front');
        end
        
    case 'equal'
        
        
        for i = 1:n_categories
            tmid = key_points(i).tmid_radians;
            
            xmid_outside = 0 + (outer_circle_radius + .04*outer_circle_radius) * cos(tmid);
            ymid_outside = 0 + (outer_circle_radius + .04*outer_circle_radius) * sin(tmid);
            
            handles(i).texth = text(xmid_outside, ymid_outside, text_labels{i}, ...
                'Rotation', myrotation(i), 'FontSize', 16 * 6/n_categories,...
                'HorizontalAlignment','center','parent',gca, 'Layer','front');
            
        end
        
    otherwise
        error('Unknown labelstyle');
        
end % labelstyle

axis off

end % main function