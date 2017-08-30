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
%        graphic that denotes a spatial brain pattern or mask created with riverplot_rect.m.  
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
%   **coveragetype:**
%        'coveragetype' followed by 'relative' 'absolute' or 'normalized'
%
%       Relative: cover 100% of each rectangle, in proportion to relative associations.
%       Thickness of ribbon is proportional to input similarity only within
%       each element, not across them
%
%       absolute coverage: Ribbon width is proportional to similarity on 0-1
%       scale, with no normalization
%
%       normalized coverage [default]: Ribbon width for layer 1 is relative to layer 1
%       so that each layer is fully covered, but on layer2 is proportional to
%       similarity across all elements, normalized by max across all
%
%       sum-normalized coverage: Ribbon width for layer 1 is relative to layer 1 
%       so that each layer is fully covered, but on layer2 is proportional to
%       summed similarity across all elements, that is normalized by max 
%       of the summ across all similarities for each layer2 item.
%
%   **color, colors:**
%        followed by cell array of {[r g b] [r g b]...} color spec for
%        ribbons. Colors are specified for each element of layer1.
%
%   **steepness**
%        followed by value for sigmoid steepness, on [0 1] inteval
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
%    8/21/2017 Stephan Geuter - added sum-normalization
% ..

% ..
%    DEFAULTS AND INPUTS
% ..

ax = gca;
[n2, n1] = size(sim_matrix);
mycolor = scn_standard_colors(n1);
coveragetype = 'normalized';
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
                
            case 'coveragetype', coveragetype = varargin{i+1}; varargin{i+1} = [];
                    
            % arguments passed on to subfunctions: Do nothing, but these
            % are used.
            case 'from_bottom'

            case 'coverage', coverage = varargin{i+1}; varargin{i+1} = [];
                
            case 'steepness', steepness = varargin{i+1}; varargin{i+1} = [];
                
            case ' ' 
                % do nothing, avoid warning from calling function
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

axes(ax)

ribbons = cell(n2, n1);

if length(layer2) ~= size(sim_matrix, 1), error('layer2 size does not match sim matrix column size'); end
if length(layer1) ~= size(sim_matrix, 2), error('layer1 size does not match sim matrix row size'); end

% For determining offset and position: 
my_total_layer1 = sum(abs(sim_matrix));
my_total_layer2 = sum(abs(sim_matrix'))';

my_total_layer1(my_total_layer1 == 0) = 1;  % prevent dividing by zero
my_total_layer2(my_total_layer2 == 0) = 1;  % prevent dividing by zero

% for both offsets and coverages. offsets are starting points, coverages
% are the thickness of each ribbon
% coverages should always be positive - take absolute values (SG 8/21/2017)
fullcoverages_layer1 = abs(bsxfun(@(x,y) x ./ y, sim_matrix, my_total_layer1));
fullcoverages_layer2 = abs(bsxfun(@(x,y) x ./ y, sim_matrix, my_total_layer2));

        
switch coveragetype
    
    case 'relative'
        % Relative: cover 100% of each rectangle, in proportion to relative associations.
        % Thickness of ribbon is proportional to input similarity only within
        % each element, not across them
        
        coverages_layer1 = fullcoverages_layer1;
        coverages_layer2 = fullcoverages_layer2;
        
    case 'absolute'
        % absolute coverage: Ribbon width is proportional to similarity on 0-1
        % scale, with no normalization
        coverages_layer1 = abs(sim_matrix);
        coverages_layer2 = abs(sim_matrix);
        
    case 'normalized'
        % normalized coverage [default]: Ribbon width for layer 1 is relative to layer 1 
        % so that each layer is fully covered, but on layer2 is proportional to
        % similarity across all elements, normalized by max across all
        
        mymax = max(abs(sim_matrix(:)));
        if mymax <=0, mymax = 1; end % in case empty
        
        %mymax = max(my_total_layer2); 
        
%         coverages_layer1 = sim_matrix ./ mymax;
        coverages_layer1 = fullcoverages_layer1;
        
        coverages_layer2 = abs(sim_matrix) ./ mymax;
        
    case 'sumnormalized'    
        % sum-normalized coverage: Ribbon width for layer 1 is relative to layer 1 
        % so that each layer is fully covered, but on layer2 is proportional to
        % summed similarity across all elements, that is normalized by max 
        % of the summ across all similarities for each layer2 item.
        
        mymax = max(sum(abs(sim_matrix)));
        if mymax <=0, mymax = 1; end % in case empty
        
        coverages_layer1 = fullcoverages_layer1;
        coverages_layer2 = abs(sim_matrix) ./ mymax;
        
    otherwise
        error('Illegal string entered for coveragetype');
        
end

% offsets are in one direction, take abs(coverages). SG 8/21/2017
offsets_layer1 = [zeros(1, n1); cumsum(abs(fullcoverages_layer1))];
offsets_layer2 = [zeros(1, n2); cumsum(abs(fullcoverages_layer2'))]';



% where there is only one target, center offsets
wh = sum(coverages_layer2 ~= 0, 2) == 1;
mycov = sum(coverages_layer2, 2) ./ 2; % half the coverage
mycov = 0.5 - mycov;  % midpoint 0.5 - half the coverage interval
offsets_layer2(wh, :) = repmat(mycov(wh), 1, size(offsets_layer2, 2));



for i = 1:n1
  
    %myoffset = i .* (1 ./ (1.5*n1));  % adjust based on number of sources to reduce overlap in ribbon endings
    
    for j = 1:n2
        
        myoffset = [offsets_layer1(j, i) offsets_layer2(j, i)]; 

        mycoverage = [coverages_layer1(j, i) coverages_layer2(j, i)];
        
        for k=1:2 %cludgey fix for riverplots that extend beyond layer 2 rectangle - this should be fixed better above
           if myoffset(k)+mycoverage(k)>1 %we dont want ribbons extending beyond range.. scoot down a bit
              myoffset(k)=max(0,myoffset(k)- (myoffset(k)+mycoverage(k)-1) ); % 8/30/17 SG scoot only down to match the max of box (=1)
           end
        end
        
        % 8/30/17 Stephan Geuter
        % reverse offsets in order to reduce the number of ribbon
        % crossings.
        myoffset = 1-myoffset-mycoverage;
        
        
        % this would be the relative coverage way
        %mycoverage = [sim_matrix(j, i) ./ my_total_layer1(i) sim_matrix(j, i) ./ my_total_layer2(j)]; 
        
%         myoffset = j .* (1 ./ (2*n2));  % adjust based on number of targets to reduce overlap in ribbon endings
        
        if sim_matrix(j, i) ~= 0
            
            ribbons{j, i} = riverplot_ribbon(layer1{i}, layer2{j}, 'coverage', mycoverage, ...
                'y_offset', myoffset, 'color', mycolor{i}, 'from_bottom', 'steepness', steepness);

        
            % riverplot_remove_ribbons(ribbons)  % for debugging.
        end
        
    end
    
end

end % function
