function myrect = riverplot_rect(varargin)
% Standard text for function documentation
%
% First line: One-line summary description of function
%
% :Usage:
% ::
%
%     [list outputs here] = function_name(list inputs here, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) <year>  <name of author>
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
%   **param1:**
%        description of param1
%
%   **param2:**
%        description of param2
%
% :Optional Inputs:
%   **param1:**
%        description of param1
%
%   **param2:**
%        description of param2
%
% :Outputs:
%
%   **out1:**
%        description of out1
%
%   **out2:**
%        description of out2
%
% :Examples:
% ::
%
% myrect = riverplot_rect('color', [.7 .3 .3], 'position', [1 2 1 3]);
%
% :References:
%   CITATION(s) HERE
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% ..
%    DEFAULTS AND INPUTS
% ..

ax = gca;
myposition = [1 3 2 5];
mycolor = [.3 .3 .8];

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'ax', ax = varargin{i+1}; varargin{i+1} = [];
            case 'position', myposition = varargin{i+1}; varargin{i+1} = [];
            case 'color', mycolor = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


recth = imrect(ax, myposition);

patchh = findobj(get(recth, 'Children'), 'Type', 'Patch');
set(patchh, 'FaceColor', mycolor);

myvert = get(patchh, 'Vertices');

bottomleft = min(myvert);  % x y coords in plot
bottomright = [max(myvert(:, 1)) min(myvert(:, 2))];
topleft = [min(myvert(:, 1)) max(myvert(:, 2))];
topright = max(myvert);

lineh = findobj(get(recth, 'Children'), 'Type', 'Line');

myrect = struct('recth', recth, 'patchh', patchh, 'lineh', lineh, 'myvert', myvert, ...
    'bottomleft', bottomleft, 'bottomright', bottomright, 'topleft', topleft, 'topright', topright);

% Turn off marker lines
set(lineh, 'Visible', 'off')

end % function



