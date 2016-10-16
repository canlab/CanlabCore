function colors = custom_colors(startrgb, endrgb, k)
% Define cell array of custom colors in a specific range of the spectrum,
% to be used in visualization.
%
% :Usage:
% ::
%     colors = custom_colors(startrgb, endrgb, k [length])
% ..
%     Author and copyright information:
%
%     Copyright (C) 2016 Tor Wager
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
%   **startrgb**
%        3-color RGB triplet, values 0-1, for starting color in range
%
%   **endrgb**
%        3-color RGB triplet, values 0-1, for starting color in range
%
%   **k**
%        how many colors you want, from 1 - 64. Colors will be sampled
%        evenly along the range you specify
%
% :Outputs:
%
%   **colors_cell**
%        Cell array of k colors
%
% :Examples:
% ::
%
%    % give examples of code here
%    param1 = abc();
%    param2 = xyz();
%    [out1,out2] = func_call(param1, param2)
%
% :See also:
%   colormap_tor, bucknerlab_colors, scn_standard_colors
%


% define 64-color map in range you specify
% ---------------------------------------------
cm = colormap_tor(startrgb, endrgb); % custom colormap

% turn it into a cell array and sample only the number desired
% ---------------------------------------------

whcolors = round(linspace(1, 64, k));

colors = mat2cell(cm(whcolors, :), ones(k, 1));

end % function


