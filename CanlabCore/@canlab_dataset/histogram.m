function fig_han = histogram(D, varargin)
% Histogram of one variable in dataset
%   - can be either event-level or subject-level
%   - event-level data is plotted as concatenated events across subject-level
%   - both variables must be valid names (case-sensitive)
%
% :Usage:
% ::
%
%    fig_han = histogram(D, [optional inputs]);
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
% :Inputs:
%
%   **D:**
%        canlab_dataset
%
% :Optional Inputs:
%   **nofig': suppress creation of new figure
%
%
% :Outputs:
%
%   **fig_han:**
%        figure handle
%


fig_han = [];
dofig = 1;

if any(strcmp(varargin, 'nofig'))
    dofig = 0;
end

[dat1, dcell1, whlevel1] = get_var(D, v1, varargin{:});


if isempty(dat1) 
    % skip
    disp('No plot: Missing variables');
    return
end

if dofig
    fig_han = create_figure([v1 ' Histogram']);
else
    fig_han = gcf;
end

% Set number of bins
nbins = max(10, length(dat1(:)) ./ 100);
nbins = min(nbins, 200);
nbins = min(nbins, length(unique(dat1(:))));

switch whlevel1
    case 1
        hist(dat1(:), nbins);
        
    case 2
        
        hist(dat1(:), nbins);
        
    otherwise
        error('Illegal level variable returned by get_var(D)');
end

        grid off
han = gca;
set(gca, 'FontSize', 24)

xlabel(v1);


end % function

