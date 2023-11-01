% Utility for table objects that extracts data from a variable, unstacks by Subject ID (or any grouping variable), and returns a cell vector and padded table with data values.
%
% :Usage:
% ::
%
%     [data_cell, data_table] = table_extract_data_by_subject(T, sid_field_name, data_field_name_to_extract)
%
% You can use this function to unstack variables in long-format tables
% (multiple rows per participant) into cell arrays and wide-format tables
% that separate observations by participant.
%
% This is useful for analyzing
% within-person relationships between variables (e.g., ratings, reaction
% times, etc. across trials within-participant), and for formatting inputs
% to be consistent with some statistical functions.  For example,
% line_plot_multisubject plots within-person relationships between
% variables and returns summary statistics like the average within-person
% correlation and stats on the relationship.
%
% Functions including mediation, mediation_brain_multilevel,
% glmfit_multilevel, and line_plot_multisubject take the output format
% returned from this function in data_cell.
%
% :Inputs:
%
%   **T:**
%        A stacked data table, with at least one variable indicating
%        subject ID (with numeric values, e.g. integers, or text)
%
%   **sid_field_name:**
%        Name of the table field indiciating subject ID
%
%   **data_field_name_to_extract:**
%        Name of the table field indiciating the data to extract
%
% :Outputs:
%
%   **data_cell:**
%        cell vector with 1 cell per subject, and a column vector of data
%        in each cell
%
%   **data_table:**
%        Table object with one column per subject, and rows indicating data
%        within-subject, padded with NaNs in case of unequal numbers of
%        observations
%
% :Examples:
% ::
%
% sid = [1 1 1 1 1 2 2 2 2 3 3 3]'; rating = [3 1 2 4 3 6 4 6 6 9 2 7]';
% T = table(sid, rating, 'VariableNames', {'sid' 'rating'});
% T
% [rating_cell, rating_table] = table_extract_data_by_subject(T, 'sid', 'rating')
%
% :See also:
%   - Use in glmfit_multilevel, multilevel_mediation, igls, and other
%   functions taking cell data as input
%

% ..
%     Author and copyright information:
%
%     Copyright (C) 2023 Tor Wager
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

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

function [data_cell, data_table] = table_extract_data_by_subject(T, sid_field_name, data_field_name_to_extract)

wasstring = ischar(T.(sid_field_name)(1));

newnames = cellstr(num2str(unique(T.(sid_field_name), 'stable')));

if ~wasstring
    newnames = cellfun(@(x) ['subj_' deblank(x)], newnames, 'UniformOutput', false);
    newnames = strrep(newnames, ' ', '');
end

% this would unstack, but problem is that (1) it aggregates, (2) may be
% unstable wrt order, (3) how to deal with diff numbers of rows is tricky
% (maybe it will?)
% subtable = T(:, {data_field_name_to_extract sid_field_name});
% U = unstack(subtable, data_field_name_to_extract, sid_field_name, 'NewDataVariableNames', newnames); %, 'GroupingVariables', data_field_name_to_extract);

% Do it manually
rownames = num2str(T.(sid_field_name));
indic = string2indicator(rownames);
indic = logical(indic);

data_cell = cell(1, length(newnames));

for i = 1:length(data_cell)

    data_cell{i} = T.(data_field_name_to_extract)(indic(:, i), :);

end

% convert to NaN-padded table
max_rows = max(sum(indic));
data_table = table;

for i = 1:length(data_cell)

    mydata = padarray(data_cell{i}, max_rows - size(data_cell{i}, 1), NaN, 'post');

    data_table.(newnames{i}) = mydata;

end

end % function






function [indic, nms, condf] = string2indicator(str,varargin)
%[indic, nms, condf] = string2indicator(str,varargin)
%
% Takes a cell vector of string labels and returns an indicator matrix
% Optional argument is a cell vector of strings for which values to match.
%
% Examples:
%[indic,nms] = string2indicator(CL{1}(1).valence);
%[indic,nms] = string2indicator(CL{1}(1).valence,{'neg' 'pos'});
%
% Note: This function changed 11/1/2023 to return names in stable order
% (the order which they appear in the original str input).
% Also, previously, the function would accept a non-cell string matrix but
% return potentially incorrect output; now enforces cell array


if ~iscell(str)
    str = cellstr(str);
end

if ~isempty(varargin)
    nms = varargin{1};
else
    nms = unique(str, 'stable');  % tor edited 11/1/2023, 'rows' would fix edge case with non-cell string matrix, but now enforce cell anyway
end

if ~iscell(nms) %, for i = 1:length(nms), nms2{i} = nms(i); end, nms = nms2; end  % this was returning only the first character of a string matrix
    nms = cellstr(nms);
end

nms(strcmp(nms, 'NaN')) = [];

indic = zeros(length(str), 1);

for k = 1:length(nms)
    indic(:,k) = strcmp(str, nms{k});
end

if size(nms,1) > size(nms,2)
    nms = nms';
end

if nargout > 2
    % make condition function
    condf = indic2condf(indic);
end

end
