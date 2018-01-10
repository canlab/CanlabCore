function fig_han = scattermatrix(D, wh_level, wh_vars)
% Scatterplot matrix of pairwise event-level variables
%
% :Usage:
% ::
%
%    fig_han = scattermatrix(D, wh_level, wh_vars)
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
%   **wh_level:**
%        1 (Subject) or 2 (Event)
%
%
% :Outputs:
%
%   **fig_han:**
%        figure handle
%
% :Examples:
% ::
%
%    fig_han = scattermatrix(D);
%
%    wh = [5:9];
%    fig_han = scattermatrix(D, 2, wh);
%
%    f = scattermatrix(D, 2, {'Choice' 'RT' 'Pain' 'SwitchNext' 'Frustration' 'Anxiety' 'Control'});
%


switch wh_level
    case 1
        names = D.Subj_Level.names;
        
    case 2
        names = D.Event_Level.names;
    otherwise
        error('Enter 1 or 2 for wh_level.');
end

if nargin < 3 || isempty(wh_vars)
    wh_vars = 1:length(names);
end

% convert from cell of names of needed
if iscell(names)
    
    myvars = get_var(D, names);
    
elseif iscell(names)
    wh = zeros(size(names));
    for i = 1:length(wh_vars)
        wh = wh + strcmp(wh_vars{i}, names);
    end
    wh_vars = find(wh);
end

names = names(wh_vars);

n = length(names);

fig_han = create_figure('scattermatrix', n, n);

for i = 1:n
    
    subplot(n, n, (i-1)*n + i)
    histogram(D, names{i}, 'nofig');
    
    for j = i+1:n
        
        subplot(n, n, (i-1)*n + j)
        
        scatterplot(D, names{i}, names{j}, 'nofig');
        
        title(' ');
        drawnow
        
    end
    
end

end % function



