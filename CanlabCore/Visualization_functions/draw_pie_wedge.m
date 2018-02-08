function [handles, key_points] = draw_pie_wedge(T1, T2, myradius, varargin)
% Draws a pie wedge with color fill, for polar plots and circular charts
%
% :Usage:
% ::
%
%     handles = draw_pie_wedge(T1, T2, myradius, ['linecolor'], [r g b], ['fillcolor'], [r g b])
%
% ..
%     Dec 2017, Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **T1:**
%        Starting location, in radians. 0 = to the right on plot; 2*pi is
%        full-circle, also to the right on plot.
%
%   **T2:**
%        Ending location, in radians. 0 = to the right on plot; 2*pi is
%        full-circle, also to the right on plot.
%
% :Optional Inputs:
%   **'linecolor':**
%        Followed by [r g b] triplet specifying line color, or 'none'
%
%   **'fillcolor:**
%        Followed by [r g b] triplet specifying fill color, or 'none'
%
%   **'outside_radius_offset'**
%        Followed by outside_radius_offset value. See below.
% 
%   **'linewidth':**
%       Followed by line width
%
%   **'stripes':**
%       Turn on radial striping; can be followed by value true or false
%
%   **'stripedensity':**
%       Number of stripes/lines to draw in wedge
%
% :Outputs:
%
%   **handles:**
%        A structure of line and fill handles
%
%   **key_points:**
%        x and y locations for the mid-point of the arc, inner mid-point of the
%        wedge, and outer mid-point offset by a constant value. This is
%        set by outside_radius_offset, and is useful for text labels
%
% :Examples:
% ::
%
% myradius = 3;
% T1 = 1; % start, in radians
% T2 = 4; % end, in radians
%
% handles = draw_pie_wedge(1, 3, 2); % Draw wedge from radian=1 to 3, radius = 2
%
% Draw reference circle with radius = 2
% handles = draw_pie_wedge(0, 2*pi, 2, 'linecolor', [0 0 0], 'fillcolor', 'none');
%
% Now draw a wedge on top of that:
% handles = draw_pie_wedge(1, 1.5, 2.5, 'linecolor', [.7 .7 0], 'fillcolor', [.3 .7 .7]);
%
% Draw with dense striping pattern:
% handles = draw_pie_wedge(1, 1.5, 2.5, 'linecolor', [.7 .7 0], 'fillcolor', [.3 .7 .7], 'stripes', 'stripedensity', 40);
%
% Draw a series of wedges in different colors:
% n_categories = 7;
% st = linspace(0, 2*pi, n_categories+1); % starting vals in radians
% 
% en = st(2:end);
% st = st(1:end-1);
% 
% colors = seaborn_colors(n_categories + 1); % add one to make colors more diff
% clear handles
% 
% figure; hold on;
% 
% for i = 1:n_categories
%     
%     handles(i) = draw_pie_wedge(st(i), en(i), i, 'linecolor', colors{i} ./ 2, 'fillcolor', colors{i});
%     
% end
% 
% axis equal
%
% For another extended example with text labels, see tor_wedge_plot.m



linecolor = [1 0 0];
fillcolor = [.7 .3 .7];
linewidth = 2;
outside_radius_offset = .04; % Value to offset by when calculating outside-wedge midpoint for labels
dostripes = false;
stripedensity = 12;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'linecolor', linecolor = varargin{i+1}; varargin{i+1} = [];
            case 'fillcolor', fillcolor = varargin{i+1}; varargin{i+1} = [];
                
            case 'linewidth', linewidth = varargin{i+1}; varargin{i+1} = [];
                    
            case 'stripes' 
                dostripes = true; 
                % Optional argument true or false following keyword
                if length(varargin) > i && ~ischar(varargin{i+1}) % param value entered in next arg
                    dostripes = varargin{i+1}; varargin{i+1} = [];
                end
                
            case 'stripedensity', stripedensity = varargin{i+1}; varargin{i+1} = [];
                
            case 'outside_radius_offset', outside_radius_offset = varargin{i+1}; varargin{i+1} = [];
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% This would draw a reference circle
%hold on
%h0 = drawCircle(0, 0, myradius);
%set(h0, 'Color', 'k');

% Circle Arc
% ------------------------------------------
%h = drawCircleArc(0, 0, myradius, T1, T2);

[h xt yt] = drawArc(T1, T2, myradius);

% get midpoint of arc
% ------------------------------------------
tmid = (T1+T2)/2; 
xmid_arc = 0 + myradius*cos(tmid);
ymid_arc = 0 + myradius*sin(tmid);

% Middle of pie
xmid_pie = 0 + myradius/2*cos(tmid);
ymid_pie = 0 + myradius/2*sin(tmid);

% outside arc
rad2 = myradius + outside_radius_offset; %.1 * myradius;
xmid_outside = 0 + rad2*cos(tmid);
ymid_outside = 0 + rad2*sin(tmid);

% others
radialangle = rad2deg(tmid); % radial angle in degrees
tangentangle = rad2deg(tmid + (pi/2)); % angle of tangent, for text

key_points = struct('tmid_radians', tmid, 'radialangle', radialangle, 'tangentangle', tangentangle, ...
    'xmid_arc', xmid_arc, 'ymid_arc', ymid_arc, 'xmid_outside', xmid_outside, ...
    'ymid_outside', ymid_outside, 'xmid_pie', xmid_pie);

set(h, 'Color', linecolor, 'LineWidth', linewidth)

% Pie wedge
% ------------------------------------------
hold on
[x1,y1] = pol2cart(T1, myradius);
h(2) = plot([0 x1], [0 y1], 'Color', linecolor, 'LineWidth', linewidth);

[x2,y2] = pol2cart(T2, myradius);
h(3) = plot([0 x2], [0 y2], 'Color', linecolor, 'LineWidth', linewidth);

% Fill
% ------------------------------------------

if ischar(fillcolor) && strcmp(fillcolor, 'none')
    % skip it
    hfill = [];
else
    hfill = fill([xt 0], [yt 0], fillcolor, 'FaceAlpha', .65);
end

handles = struct('line_han', h, 'fill_han', hfill);

% Striping
% ------------------------------------------
if dostripes
    
    handles.stripe_han = drawStripes(T1, T2, myradius, stripedensity);
    set(handles.stripe_han, 'Color', linecolor, 'LineWidth', linewidth)
    
else
    
    handles.stripe_han = [];
    
end
    


end % main function


function [h xt yt] = drawArc(T1, T2, myradius)

% key bits
t = T1:.01:T2;          % Arc
xt = 0 + myradius*cos(t);
yt = 0 + myradius*sin(t);
xt = [xt 0+myradius*cos(T2)];
yt = [yt 0+myradius*sin(T2)];

h = line(xt, yt);

end


function h = drawStripes(T1, T2, myradius, stripedensity)

r = linspace(0, myradius, stripedensity);

for i = 1:length(r)
    
    h(i) = drawArc(T1, T2, r(i));
    
end

end

