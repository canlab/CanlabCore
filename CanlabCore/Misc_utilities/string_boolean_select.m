function [cell_indices, uniquelabels] = string_boolean_select(input_cell_str, strings_to_find, varargin)
% Select indices from a cell array of strings that match any of a series of input cell string stems
%
% [cell_indices, uniquelabels] = string_boolean_select(input_cell_str, ['not', strings_to_exclude], ['and', strings_to_conjoin])
%
% - Select indices with text match any pattern (char string) from a cell array of strings (input_cell_str)
% - If 'not', next input is a cell array of strings defining patterns to exclude (strings_to_exclude)
% - If 'and', next input is a cell array of strings defining patterns to conjoin with and (strings_to_conjoin)
%   First, any cells matching and of the stems in strings_to_conjoin are
%   identified, then logical and is performed.
%
% :Inputs:
%
%   **input_cell_str:**
%        Cell vector of strings with text to match
%
%   **strings_to_find:**
%        Cell vector of strings with patterns to match,
%        e.g., {'HR', 'HF'} would find any strings in input_cell_str that
%        match either, e.g., 'HRV' and 'MYHF' would be a hits.
%
% :Optional Inputs:
%   **'not':**
%        Followed by cell vector of strings with patterns to exclude (strings_to_exclude).
%        Any cell in strings_to_find matching any of the patterns in strings_to_exclude
%        is excluded. e.g., if strings_to_find is {'HR', 'HF'} is
%        strings_to_exclude is {'F'}, 'HRV' would match but 'MYHF' would
%        not.
%
%   **'and':**
%        Followed by cell vector of strings with patterns to conjoin (strings_to_conjoin)
%        e.g., if strings_to_find is {'HR', 'HF'} and strings_to_conjoin is
%        {'F'}, 'MYHF' would match but 'HRV' would not.
%
% :Outputs:
%
%   **cell_indices:**
%        Vector if which cells in input_cell_str match
%
%   **uniquelabels:**
%        Unique values of matching cells in input_cell_str
%

% to_exclude = false(size(input_cell_str));  % 'not'
% to_conjoin = false(size(input_cell_str));  % 'and'

strings_to_exclude = {};
strings_to_conjoin = {};


% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'not', strings_to_exclude = varargin{i+1}; varargin{i+1} = [];
            case 'and', strings_to_conjoin = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% -------------------------------------------------------------------------
% FIND BY STRING
% -------------------------------------------------------------------------
cell_indices = match_any_pattern(input_cell_str, strings_to_find);


% NOT
if ~isempty(strings_to_exclude)
    
    to_exclude = match_any_pattern(input_cell_str, strings_to_exclude);
    cell_indices(to_exclude) = false;
    
end

% AND

if ~isempty(strings_to_conjoin)
    
    to_conjoin = match_any_pattern(input_cell_str, strings_to_conjoin);
    cell_indices = cell_indices & to_conjoin;
    
end

uniquelabels = unique(input_cell_str(cell_indices));

end % main function





function cell_indices = match_any_pattern(input_cell_str, strings_to_find)

cell_indices = false(size(input_cell_str));  % main search

for i = 1:length(strings_to_find)
    
    % Find which names match
    wh = ~cellfun(@isempty, strfind(input_cell_str, strings_to_find{i}));
    
    cell_indices = cell_indices | wh;
    
end

end
