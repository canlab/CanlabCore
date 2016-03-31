function [h, p, ci, stats] = ttest2(D, varname, wh_keep1, wh_keep2, varargin)
% Two sample ttest for two samples of one subject-level variable
%
% :Usage:
% ::
%
%    ttest2(D, varname, wh_keep1, wh_keep2, [optional inputs])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2013 Tor Wager
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
%
% :Inputs:
%
%   **D:**
%        a canlab_dataset object
%
%   **varname:**
%        the name of a valid variable to get from dataset
%
%   **wh_keep1:**
%        subjects forming first sample              
%
%   **wh_keep2:**
%        subjects forming second sample              
%
% :Optional Inputs:
%
%   **noverbose:**
%         will suppress print out of results and bargraph
%
%   **varargin:**
%         other variables passed directly to MATLAB's ttest2
%
% :Outputs:
%
%   same as MATLAB's ttest2 output 
%


if any(wh_keep1 & wh_keep2), warning('YOUR SAMPLES ARE OVERLAPPING!!'); end

verbose=1;
if any(strcmp('noverbose', varargin))
    verbose=0;
    varargin(find(strcmp('noverbose', varargin))) = [];
end

x1 = get_var(D, varname, wh_keep1);
x2 = get_var(D, varname, wh_keep2);

[h, p, ci, stats] = ttest2(x1, x2, varargin{:});

if verbose, ttest2_printout(x1,x2, 1); end

end
