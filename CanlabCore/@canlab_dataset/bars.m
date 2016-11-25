function [dat, descrip, colors, h1, s1] = bars(obj, varnames, varargin)
% Bar plot for canlab_dataset object
%
% :Usage:
% ::
%
%    [dat, descrip, colors, h1, s1] = bars(obj, varnames, [optional inputs])
%
%    Takes inputs to barplot_columns
%    e.g., 'noviolin' to suppress violin plot
%    'colors' followed by colors cell for each bar
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
%   **obj:**
%        canlab_dataset object
%
%   **varnames:**
%        Cell string of variable names to plot
%
% :Optional Inputs:
%    **colors:**
%        defined colors (default are set by scn_standard_colors.m)
%
%   **nofig:**
%        do not generate figure
%
%   Takes any optional inputs to barplot_colored.m
%
%
% :Outputs:
%
%   **dat:**
%        data matrix for each variable
%
%   **descrip:**
%        the description for this variable
%
%   **colors:**
%        selected colors (default are set by scn_standard_colors.m)
%
%   **h1:**
%        figure handle
%
%   **s1:**
%        axis handle
%
% :Examples:
% ::
%
%     create_figure('NPS values - All subjects');
%
%     varnames = {'15' '13' '11' ' 9' '16' '14' '12' '10'};
%     xvals = [1 2 4 5 8 9 11 12];
%     colors = {[1 0 0] [0 1 0] [1 0 0] [0 1 0] [1 0 0] [0 1 0] [1 0 0] [0 1 0]};
%     bars(LevoNPS, varnames, 'x', xvals, 'colors', colors, 'XTickLabels', varnames, 'within', 'nofig');
%

n = length(varnames);
colors = scn_standard_colors(n);


% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% optional inputs with default values
% -----------------------------------
% - allowable_args is a cell array of argument names
% - avoid spaces, special characters, and names of existing functions
% - variables will be assigned based on these names
%   i.e., if you use an arg named 'cl', a variable called cl will be
%   created in the workspace

allowable_args = {'colors' 'nofig'};

default_values = {colors [0]};

% define actions for each input
% -----------------------------------
% - cell array with one cell for each allowable argument
% - these have special meanings in the code below
% - allowable actions for inputs in the code below are: 'assign_next_input' or 'flag_on'

actions = {'assign_next_input' 'flag_on'};

% logical vector and indices of which inputs are text
textargs = cellfun(@ischar, varargin);
whtextargs = find(textargs);

for i = 1:length(allowable_args)
    
    % assign default
    % -------------------------------------------------------------------------
    
    eval([allowable_args{i} ' =  default_values{i};']);
    
    wh = strcmp(allowable_args{i}, varargin(textargs));
    
    if any(wh)
        % Optional argument has been entered
        % -------------------------------------------------------------------------
        
        wh = whtextargs(wh);
        if length(wh) > 1, warning(['input ' allowable_args{i} ' is duplicated.']); end
        
        switch actions{i}
            case 'assign_next_input'
                eval([allowable_args{i} ' = varargin{wh(1) + 1};']);
                %varargin{wh(1) + 1} = [];
                
            case 'flag_on'
                eval([allowable_args{i} ' = 1;']);
                
            otherwise
                error(['Coding bug: Illegal action for argument ' allowable_args{i}])
        end
        
    end % argument is input
end

% END DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% REST OF BAR PLOT


[dat, ~, ~, descrip] = get_var(obj, varnames, varargin{:});

%descrip = cat(1, descrip{:});
%colors = cat(1, colors{:});

xvals = 1:n;

if nofig
    % Do nothing, use existing figure
else
    create_figure('bars');
end

[h1, s1] = barplot_columns(dat, 'nofig', varargin{:});
% set(h2, 'BarWidth', .9)
%colormap(colors)

mynames = strrep(varnames, '_', ' ');
set(gca, 'XTickLabel', mynames);

% STATS OUTPUT

fprintf('T-tests against zero\n');
fprintf('Description\tt(%3.0f)\tp\n', size(dat, 1) - 1);

% Print table - t-test against zero
for i = 1:n
    [h, p, ci, stats] = ttest(dat(:, i));
    if p < .001, str = '***';
    elseif p < .01, str = '**';
    elseif p < .05, str = '*';
    elseif p < .10, str = '+';
    else str = ' ';
    end
    
    fprintf('%s\t%3.2f\t%3.6f\t%s\n', descrip{i}, stats.tstat, p, str);
end


end

