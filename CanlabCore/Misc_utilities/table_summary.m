function summaryTable = table_summary(table_obj)
% table_summary Create a summary table of variable characteristics.
%
% :Usage:
% ::
%     summaryTable = table_summary(table_obj)
%
% :Inputs:
%
%   **table_obj:** 
%        A table object containing one or more variables.
%
% :Outputs:
%
%   **summaryTable:** 
%        A table with one row per variable of table_obj. The summary includes:
%
%        - **Variable:**  
%             The name of the variable.
%
%        - **DataType:**  
%             The class (data type) of the variable.
%
%        - **UniqueCount:**  
%             The number of unique values in the variable, if it is numeric, 
%             a character array, a string array, or a cell array of char. 
%             Otherwise, NaN.
%
%        - **missingCases:**  
%             The number of missing cases in the variable (computed using ismissing).
%
% :Examples:
% ::
%     % Example: Load a table from a CSV file and create a summary.
%     T = readtable('mydata.csv');
%     summaryTable = table_summary(T);
%     disp(summaryTable);
%
% :See also:
%   ismissing, unique, readtable, writetable
%
% Author: Your Name
% Date: YYYY-MM-DD
% License: GNU General Public License v3 or later

    % Get the variable names and determine the number of variables
    varNames = table_obj.Properties.VariableNames;
    nVars = numel(varNames);
    
    % Preallocate cell arrays and vectors for the summary information.
    dataTypes = cell(nVars, 1);
    uniqueCounts = nan(nVars, 1);
    missingCases = nan(nVars, 1);
    
    % Loop over each variable in the table.
    for i = 1:nVars
        varData = table_obj.(varNames{i});
        
        % Determine the data type.
        dataTypes{i} = class(varData);
        
        % Count missing cases (using ismissing works for many data types)
        missingCases(i) = sum(ismissing(varData));
        
        % Count unique values.
        if isnumeric(varData)
            uniqueCounts(i) = numel(unique(varData));
        elseif ischar(varData)
            uniqueCounts(i) = numel(unique(cellstr(varData)));
        elseif isstring(varData)
            uniqueCounts(i) = numel(unique(varData));
        elseif iscell(varData) && all(cellfun(@ischar, varData))
            uniqueCounts(i) = numel(unique(varData));
        else
            uniqueCounts(i) = NaN;
        end
    end
    
    % Create a summary table with the collected information.
    summaryTable = table(varNames', dataTypes, uniqueCounts, missingCases, ...
        'VariableNames', {'Variable', 'DataType', 'UniqueCount', 'missingCases'});
end
