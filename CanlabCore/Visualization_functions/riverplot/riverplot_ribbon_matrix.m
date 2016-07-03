function ribbons = riverplot_ribbon_matrix(layer1, layer2, sim_matrix, varargin)
%
% This is a subfunction of the riverplot toolbox.
%  - Draws a set of ribbon FROM layer1 TO layer2
%  - layers are structures created by riverplot_draw_layer.m
%  - Ribbon covers a specified percentage of each rect structure
%  (rectangle) within each layer, can be determined by a symmetrical
%  similarity matrix or an asymmetric one (sim_matrix).
%
% :Usage:
% ::
%
%     myribbon = riverplot_ribbon_matrix(layer1, layer2, sim_matrix, [optional inputs])
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
%   **layer1:**
%        Cell array of structures with handles and decription of imrect class objects,
%        created with riverplot_draw_layer.m.  Each object is a rectangular
%        graphic that denotes a spatial brain pattern or mask created with riverplot_rect.m.  This is always the LEFTMOST
%        layer1 is the LEFTMOST layer in the plot.
%
%   **layer2:**
%        Same as layer1, but this is the RIGHTMOST layer of the pair
%
% 
%   **sim_matrix:**
%       A [j x k] rectangular matrix of similarities to plot, where 
%       j is num. of elements in layer2, and k is num. elements in layer1
%       Any similarity metric of choice is OK, but all values should be between
%       0 and 1], so negative values cannot be used.
%
%
% :Optional Inputs:
%   **coverage:**
%        'coverage' followed by 2-element vector for [rect1 rect2], where
%        each value is between 0 and 1 and specifies the percentage of the
%        plot the ribbon will cover.
%
%   **color, colors:**
%        followed by cell array of {[r g b] [r g b]...} color spec for
%        ribbons. Colors are specified for each element of layer1.
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
[n2, n1] = size(sim_matrix);
mycolor = scn_standard_colors(n1);
layer1fullcoverage = false;
steepness = .05;

% arguments passed on to subfunctions:
% from_bottom

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'ax', ax = varargin{i+1}; varargin{i+1} = [];
            %case 'position', myposition = varargin{i+1}; varargin{i+1} = [];
            case {'color', 'colors'}, mycolor = varargin{i+1}; varargin{i+1} = [];
                
            % arguments passed on to subfunctions: Do nothing, but these
            % are used.
            case 'from_bottom'

            case 'layer1fullcoverage', layer1fullcoverage = true;
                
            case 'steepness', steepness = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

axes(ax)

ribbons = cell(n2, n1);

for i = 1:n1
  
    myoffset = i .* (1 ./ (1.5*n1));  % adjust based on number of sources to reduce overlap in ribbon endings

    for j = 1:n2
        
        if layer1fullcoverage
            
            mycoverage = [1 sim_matrix(j, i)];
        else
            
            mycoverage = [sim_matrix(j, i) sim_matrix(j, i)];   
        end
        
%         myoffset = j .* (1 ./ (2*n2));  % adjust based on number of targets to reduce overlap in ribbon endings
        
        if sim_matrix(j, i) ~= 0
            
            ribbons{j, i} = riverplot_ribbon(layer1{i}, layer2{j}, 'coverage', mycoverage, ...
                'y_offset', myoffset, 'color', mycolor{i}, 'from_bottom', 'steepness', steepness);
        
        end
        
    end
    
end

end % function
