function layerstruct = riverplot_layer_names(layerstruct, layernames, varargin)
%
% This is a subfunction of the riverplot toolbox.
%  - Adds or removes names from a layer struct, and corresponding text on plot
%  - Toggles between adding and removing, depending on whether names exist
%  in layerstruct.  Layerstruct is created with riverplot_draw_layer.
%
% :Usage:
% ::
%
%     layerstruct = riverplot_layer_names(layerstruct, layernames, [optional inputs])
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
%   **layerstruct:**
%        Structure with handles and decription of imrect class objects for
%        a layer. Layerstruct is created with riverplot_draw_layer.
%
%   **layernames:**
%        A cell array of strings for names
%
% :Optional Inputs:
%
%
% :Outputs:
%
%   **layerstruct:**
%        structure with names and name info added
%
%
% :Examples:
% ::
%
%
% :References:
%   
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    Created by Tor Wager, July 2, 2016
%
%    8/21/2017 Stephan Geuter
%       - updated text alignment
% ..

% ..
%    DEFAULTS AND INPUTS
% ..

ax = gca;
colors = custom_colors([.2 .2 .2], [.2 .2 .2], length(layerstruct)); % series of cells with colors for each layer
leftright_string = 'left';

if ~iscell(layernames), layernames = cellstr(layernames); end

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'ax', ax = varargin{i+1}; varargin{i+1} = [];
            %case 'position', myposition = varargin{i+1}; varargin{i+1} = [];
            case {'color', 'colors'}, colors = varargin{i+1}; varargin{i+1} = [];
                
            case 'left', leftright_string = 'left';

            case 'right', leftright_string = 'right';

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

axes(ax)

n_rects = length(layerstruct);

if n_rects > length(layernames)
    error('Not enough names entered to riverplot_layer_names.m');
end

for i = 1:n_rects
    
    % Get y position: midpoint of layer
    yval = layerstruct{i}.topleft(2) - .5 * (layerstruct{i}.topleft(2) - layerstruct{i}.bottomleft(2));
    
    % Get x position: 
    switch leftright_string
        case 'left'
            
            xval = layerstruct{i}.topleft(1) - .3;
    
            hh = text(xval, yval, layernames{i}, 'Color', colors{i}, 'FontSize', 16,...
                      'HorizontalAlignment','right');
            
        case 'right'
       
            xval = layerstruct{i}.topright(1) + .3;
            
            hh = text(xval, yval, layernames{i}, 'Color', colors{i}, 'FontSize', 16,...
                      'HorizontalAlignment','left');
    end

% removed because text property 'HorizontalAlignment' is better and adjusts 
% after re-sizing the figure. SG 8/21/2017    
%     hh = text(xval, yval, layernames{i}, 'Color', colors{i}, 'FontSize', 16);
%     
%     switch leftright_string
%         case 'left'
%             hh = shift_text(hh);
%     end
    
    layerstruct{i}.name_handles(i) = hh;
    layerstruct{i}.names(i) = layernames(i);
    
end


end % function



% function text_handle = shift_text(text_handle)
% % Center text - Adjust for width of text string
% % Given text handle, subtract 1/2 the width of text string from x and y positions
% % Borrowed/adapted from John Barber, http://www.mathworks.com/matlabcentral/fileexchange/30671-calcticks
% 
% textExt = get(text_handle, 'Extent');
% textHeight = textExt(4);
% textWidth = textExt(3);
% 
% % If using a proportional font, shrink text width by a fudge factor to
% % account for kerning.
% % ax = gca;
% % 
% % if ~strcmpi(ax.FontName,'FixedWidth')
% %     textWidth = textWidth*0.8;
% % end
% 
% mypos = get(text_handle, 'Position');
% mypos(1) = mypos(1) - textWidth;
% set(text_handle, 'Position', mypos);
% 
% end
