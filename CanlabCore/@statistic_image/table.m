function r = table(t, varargin)
%
% Create and print a table of labeled regions from thresholded statistic_image object t.
%
% :Usage:
% ::
%
%     r = table(t, [optional inputs])
%
% - Returns r, a region class object with one element per activation blob
% - r(i).shorttitle contains a label for blob i, with labels assigned according to an atlas class object (using a default atlas)
% - r(i).title contains more information about coverage of multiple labeled atlas regions.
% - Passes through optional inputs to region.table
%
% ..
%
% :Inputs:
%
%   **t:**
%        A statistic_image class object
%
% :Optional Inputs:
%  
%  See region.table for any possible pass-through optional inputs
%
% :Outputs:
%
%   **r:**
%        a region class object with one element per activation blob
%        r(i).shorttitle contains a label for blob i, with labels assigned according to an atlas class object (using a default atlas)
%        r(i).title contains more information about coverage of multiple labeled atlas regions.
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
%   region.table
%   get_wh_image, to select a single image if t contains multiple images
%


% ..
%     Author and copyright information:
%
%     Copyright (C) 2021 Tor Wager
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
%    Programmers' notes:
%   2021-3-9 Created by Tor Wager 
% ..

if size(t.dat, 2) > 1
    
    disp('table() can only handle objects with a single image. use t = get_wh_image(t, i) to select image i first.');
    
    error('Exiting');
    
end

[r_pos, r_neg] = table(region(t));

r = [r_pos r_neg];

end
