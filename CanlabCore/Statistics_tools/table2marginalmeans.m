function [margmeans, out_dat, levelnames] = table2marginalmeans(data_table, var_to_extract, varargin)
% Extract marginal means and cell array of data for a combination of levels in a table() object
%
% :Usage:
% ::
%
%     [margmeans, out_dat, levelnames] = table2marginalmeans(data_table, var_to_extract, [optional inputs])
%
% :Inputs:
%
%   **data_table:**
%        Table object with variables named in data_table.Properties.VariableNames
%
%   **var_to_extract:**
%        name of a variable in the dataset to extract means for
%
% :Optional Inputs:
%   **'vars':**
%        Followed by a cell array of variable names to condition on
%        Combinations of the conditioning variables will be used to extract
%        data and mean values for each unique combination.
%        The intention is that these should be categorical variables with
%        only a few unique values.
%        e.g., if 'drug' has 2 levels and 'regul' has 3 levels, data will
%        be extracted for 6 cells, i.e., the 2 x 3 combination.
%
%   **'exclude':**
%        followed by a logical vector of rows to exclude.
%        Num rows must match table rows.
%
% :Outputs:
%
%   **margmeans:**
%        Marginal means
%
%   **out_dat:**
%        Cell array of observations in each cell
%        Suitable for plotting individual observations or further analysis
%
%   **levelnames:**
%        Names of levels of each combination
%
% :Examples:
% ::
%
% [margmeans, out_dat, levelnames] = table2marginalmeans(data_table, 'pexp_val', 'vars', {'drug' 'regul'})
% barplot_columns(out_dat, 'names', levelnames)
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%     Author and copyright information:
%
%     Copyright (C) 2020 Tor Wager
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


% Defaults

vars = data_table.Properties.VariableNames;
exclude = [];

allowable_inputs = {'vars' 'exclude'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            
            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(var_to_extract)
    error('Enter the name of a variable in the dataset to extract means for.')
end

% Get the dataset with variables to match on

data_to_match = [];
allnames = data_table.Properties.VariableNames;
wh_vars = false(1, length(allnames));

for i = 1:length(vars)
    wh_vars(strcmp(allnames, vars{i})) = true;
end

data_to_match = data_table(:, wh_vars);

if ~isempty(exclude)
    data_to_match = data_to_match(~exclude, :);
end

% redefine vars to match columns
vars = allnames(wh_vars);

% Get unique combinations and extract data for each

urows = unique(data_to_match, 'rows');

for i = 1:size(urows, 1)
    
    wh = ismember(data_to_match, urows(i, :), 'rows'); % rows unnecessary, calls table.ismember
    
    out_dat{i} = data_table.pexp_val(wh);
    
    % Name this level
    levelnames{i} = '';
    for j = 1:length(vars)
        if isnumeric(data_table.(vars{j}))
            levelnames{i} = [levelnames{i} sprintf('%s=%d ', vars{j}, table2array(urows(i, j)))];
        else
            % Needs more work for char array variables!
        end
    end
    
end

margmeans = cellfun(@nanmean, out_dat);

end % main function