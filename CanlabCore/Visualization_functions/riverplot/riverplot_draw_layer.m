function mylayer = riverplot_draw_layer(x_loc, n_rects, varargin)
%
% This is a subfunction of the riverplot toolbox.
%  - Draws a series of rects at an x location x_loc, number specified by n_rects
%  - rect1 and rect2 are structures created by riverplot_rect.m
%  - Ribbon covers a specified percentage of each rect
%  - Can start drawing from the top or bottom of each rect
%
% :Usage:
% ::
%
%     myribbon = riverplot_ribbon(rect1, rect2, [optional inputs])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2016  Tor Wager
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
%   **rect1:**
%        Structure with handles and decription of imrect class object,
%        created with riverplot_rect.m.  This is always the LEFTMOST
%        rect of the pair.
%
%   **rect2:**
%        Same as rect1, but this is the RIGHTMOST rect of the pair
%
% :Optional Inputs:
%   **coverage:**
%        'coverage' followed by 2-element vector for [rect1 rect2], where
%        each value is between 0 and 1 and specifies the percentage of the
%        plot the ribbon will cover.
%
%   **color:**
%        followed by [r g b] color spec for ribbon
%
%   **ax:**
%        followed by axis handle to plot
%
%   **from_bottom:**
%        align ribbon with bottom edges of rects 
%
% :Outputs:
%
%   **myribbon:**
%        structure with handles and info for ribbon
%
%
% :Examples:
% ::
%
% create_figure('riverplot'); set(gca, 'YDir', 'reverse');
% layer1 = riverplot_draw_layer(2, 5);
% layer2 = riverplot_draw_layer(5, 3, 'colors', bucknerlab_colors(3));
% ribbon1 = riverplot_ribbon(layer1{1}, layer2{2}, 'coverage', [.5 .2], 'color', [.2 .5 .2], 'from_bottom');
% ribbon2 = riverplot_ribbon(layer1{1}, layer2{3}, 'coverage', [1 .5], 'color', [.7 .5 .5]);
%
% :References:
%   CITATION(s) HERE
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    Created by Tor Wager, July 2, 2016
% ..

% ..
%    DEFAULTS AND INPUTS
% ..

ax = gca;
colors = custom_colors([.4 .2 0], [.8 .5 .2], n_rects); % series of cells with colors for each layer

y_loc = 0;
x_wid = .5;
y_wid = 1;


% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'ax', ax = varargin{i+1}; varargin{i+1} = [];

            case {'color', 'colors'}, colors = varargin{i+1}; varargin{i+1} = [];
                
            case 'from_bottom', from_bottom = true;

            case 'coverage', coverage = varargin{i+1}; varargin{i+1} = [];

            case {'y_loc', 'y_start'}, y_loc = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

axes(ax)

% set y starting point for each
totalwid = y_wid + .1*y_wid;     % y space occupied, total
y_loc = y_loc:totalwid:(y_loc + totalwid.* n_rects);

if length(colors) < n_rects
    error('Too few colors entered!');
end

for i = 1:n_rects
    mylayer{i} = riverplot_rect('color', colors{i}, 'position', [x_loc y_loc(i) x_wid y_wid]);
end


end % function