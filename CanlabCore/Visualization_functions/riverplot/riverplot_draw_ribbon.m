function myribbon = riverplot_draw_ribbon(topleft, topright, bottomleft, bottomright, varargin)

%
% This is a subfunction of the riverplot toolbox.
%  - Draws a ribbon given four [x y] coordinate pairs specifying the
%  corners of the ribbon in clockwise order.  
%  - Used in riverplot_ribbon.m
%
% :Usage:
% ::
%
%     myribbon = riverplot_draw_ribbon((topleft, topright, bottomleft, bottomright, [optional inputs])
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
%   **steepness**
%        followed by value for sigmoid steepness, on [0 1] inteval
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
% rect1 = riverplot_rect('color', [.7 .3 .3], 'position', [1 2 1 3]);
% rect2 = riverplot_rect('color', [.3 .3 .6], 'position', [5 4 1 3]);
% newbr = rect2.topleft;
% newbr(2) = newbr(2) - 1.5;
% myribbon = riverplot_draw_ribbon(rect1.topright, rect2.topleft, rect1.bottomright, newbr, 'color', [.7 .3 .3]);
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
%myposition = [1 3 2 5];
mycolor = [.3 .3 .8];
steepness = .05;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'ax', ax = varargin{i+1}; varargin{i+1} = [];
            %case 'position', myposition = varargin{i+1}; varargin{i+1} = [];
            case 'color', mycolor = varargin{i+1}; varargin{i+1} = [];
                
            case 'steepness', steepness = varargin{i+1}; varargin{i+1} = [];
                
            case {'coverage', 'from_bottom', 'y_offset'} % pass it on, irrelevant but passed from higher level calling fcn
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

axes(ax)

% each line coordinate is an [x,y] pair, e.g., topleft = [x,y] for top left
% point in absolute plot coordinates

line1 = riverplot_line(topleft, topright, 'w', 1, steepness);

line2 = riverplot_line(bottomleft, bottomright, 'w', 1, steepness);

ph = patch([line1.xcoords; line2.xcoords(end:-1:1)], [line1.ycoords; line2.ycoords(end:-1:1)], mycolor);

myribbon = struct('line1', line1, 'line2', line2, 'patchh', ph);

set(line1.h, 'Color', mycolor);
set(line2.h, 'Color', mycolor);
set(ph, 'FaceAlpha', .4);


end % function