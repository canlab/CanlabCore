function myribbon = riverplot_ribbon(rect1, rect2, varargin)
%
% This is a subfunction of the riverplot toolbox.
%  - Draws a ribbon FROM rect1 TO rect2
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
%   **y_offset:**
%        'y_offset' followed by one or 2-element vector for [offset1 offset2], where
%        each value is between 0 and 1 and specifies the percentage of the
%        square for layer 1 and layer 2 to shift down before plotting ribbon.
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
% myribbon = riverplot_ribbon(rect1, rect2);
% myribbon = riverplot_ribbon(rect1, rect2, 'coverage', [1 .5], 'color', [.2 .5 .2]);
% myribbon = riverplot_ribbon(rect1, rect2, 'coverage', [.5 .2], 'color', [.2 .5 .2], 'from_bottom');
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
from_bottom = false;
coverage = [1 1];  % percent of rect 1 and rect 2 that are covered.

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'ax', ax = varargin{i+1}; varargin{i+1} = [];
                
            case 'from_bottom', from_bottom = true;
                
            case 'coverage', coverage = varargin{i+1}; varargin{i+1} = [];
                
            case 'y_offset', y_offset = varargin{i+1}; varargin{i+1} = [];
                
            % Just pass these vars on
                
            case 'position' %myposition = varargin{i+1}; varargin{i+1} = [];
            case {'color', 'colors'} %mycolor = varargin{i+1}; %varargin{i+1} = [];
            case 'steepness' % just pass it on
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

axes(ax)

% left is rect 1, right is rect 2
tl = rect1.topright;    % start point for line 1 - this is the x,y pair for left-hand rectangle, top right corner
tr = rect2.topleft;     % start point for line 1 - this is the x,y pair for right-hand rectangle
bl = rect1.bottomright; % start point for line 2
br = rect2.bottomleft;  % start point for line 2

ylen1 = rect1.topright(2) - rect1.bottomright(2);  % total y distance, left rect
ylen2 = rect2.topright(2) - rect2.bottomright(2);  % total y distance, right rect

% y offset: shift down if specified and we don't have full coverage
% based on percentage of total rect y length
% If we are covering the whole interval then ignore.
if length(y_offset) == 1, y_offset = [y_offset y_offset]; end

if coverage(1) < 1, offset1 = y_offset(1) .* ylen1; else offset1 = 0; end
if coverage(2) < 1, offset2 = y_offset(2) .* ylen2; else offset2 = 0; end

% adjust distance for coverage
ylen1 = coverage(1) * ylen1;
ylen2 = coverage(2) * ylen2;

if from_bottom
    
    % adjust start point
     % tl(1) and tr(1) are ok, they define x for first line
    % tl(2) and tr(2) need y coords adjusted by offset values 
    if any(y_offset) 
        
        tl(2) = tl(2) - offset1;
        
        tr(2) = tr(2) - offset2;
        
    end
    
    % Adjust end point based on coverage
    % tl and tr are ok, they define first line
    % bl and br need to have y (2) adjusted,but not x (1)
    % rect 1 adjust y
    bl(2) = tl(2) - ylen1; % start with tl, tr line. adjust y coord for bl relative to tl
    
    % rect 2 adjust y
    br(2) = tr(2) - ylen2; % start with tl, tr line. adjust y coord for br relative to tr
    
else
    
    % adjust start point
    if y_offset 
        
        tl(2) = tl(2) - offset1;
        
        tr(2) = br(2) - offset2;
        
    end
    
    % Adjust end point based on coverage
    
    % rect 1 adjust y
    bl(2) = tl(2) - ylen1;
    
    % rect 2 adjust y
    br(2) = tr(2) - ylen2;
    


end

% will draw two lines: from [x,y] pair tl to [x,y] pair tr, and same for bl, br
myribbon = riverplot_draw_ribbon(tl, tr, bl, br, varargin{:});


end % function