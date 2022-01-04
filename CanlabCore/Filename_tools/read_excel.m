function DAT = read_excel(excelfilename, has_header_row)
% Reads each sheet of an Excel input file into a data structure DAT
%
% :Usage:
% ::
%
%     DAT = read_excel(excelfilename, has_header_row [1 or 0])
%
%  - reads multiple sheets
%  - adds fields named with variable names for easy access
%  - uses importdata <matlab internal> to do most of the work
%
% ..
%     Author and copyright information:
%     -------------------------------------------------------------------------
%     Copyright (C) 2015 Tor Wager
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
% :Input:
%
%   **excelfilename:**
%        Char array with full path of input filename
%
% :Outputs:
%
%   **DAT:**
%        Structure with DAT.(sheetname).(varname) data fields
%        DAT.(sheetname).varnames is cell array of var names
%
%   **behdat:**
%        Data structure as read in by importdata
%
% :Examples:
% ::
%
%    excelfilename = fullfile(basedir, 'data', 'SchulzNewcorn_regressor_demograph.xlsx');
%    has_header_row = 1;
%    [DAT, behdat] = read_excel(excelfilename, has_header_row)
%
% :See also:
% importdata, read_database, read_database2, read_edat_output_2008, read_physio_data
%
% ..
%    Programmers' notes:
%    Created by Tor, 7/2015
%    Updated by Tor, Dec 2021 - keep up with changes in importdata behavior
% ..

DAT = struct;

% Note: importdata is not handling blanks in Excel files well. Use xlsread below.
behdat = importdata(excelfilename); % read the original data

has_multiple_sheets = isstruct(behdat.data) || isstruct(behdat.textdata);

if has_multiple_sheets

    % % Sheet names
    % %-------------------
    mysheets = fieldnames(behdat.data);

    for s = 1:length(mysheets)

        [~, ~, raw] = xlsread(excelfilename, s);  % 2021: it doesn't extract data well with blank entries. Use raw and do it ourselves.

        [data_struct, data_table] = process_sheet(raw, has_header_row);
        
        DAT.(mysheets{s}).raw = raw;
        DAT.(mysheets{s}).data_struct = data_struct;
        DAT.(mysheets{s}).data_table = data_table;

        % DAT.(mysheets{s}).data_table = process_sheet(behdat.data.(mysheets{s}), behdat.textdata.(mysheets{s}), has_header_row);

    end

else
    % One sheet
    [~, ~, raw] = xlsread(excelfilename, 1);

    DAT.raw = raw;
    [DAT.data_struct, DAT.data_table] = process_sheet(raw, has_header_row);

    %     DAT.data_table = process_sheet(behdat.data, behdat.textdata, has_header_row);
end

end % main function


% ------------------------------------------------------------------------
% Subfunctions
% ------------------------------------------------------------------------


function [data_struct, data_table] = process_sheet(raw, has_header_row)
% Create a table object from raw data from xlsread

% % Variable names
if has_header_row

    varnames = raw(1, :);

else
    varnames = cell(1, size(raw, 2));  % textdata only seems to be full-width

    for i = 1:length(varnames)
        varnames{i} = sprintf('Var_%03d', i);
    end

end

varnames = matlab.lang.makeValidName(varnames);

% % Loop through columns
% --------------------

data_struct = struct();

for i = 1:length(varnames)

    fname = varnames{i};

    % figure out if is text and process (replace empty)
    mydat = raw((has_header_row+1):end, i); % skip first row if headers.

    % replace empty, convert to numeric if possible
    mydat = process_data_column(mydat);


    data_struct.(fname) = mydat;

end

data_table = struct2table(data_struct);

end



function mydat = process_data_column(mydat)

% Are some characters text?
is_char = any(cellfun(@ischar, mydat));

if ~is_char
    % already numeric. convert from cell.
    mydat = cat(1, mydat{:});
    return
end

% Below, some entries at least are char arrays. But may be numeric with
% missing values/empty cells

% some empty cells in text are read in as NaN, which can cause errors below
nanvec = false(size(mydat));
for i = 1:length(mydat)
    if isnan(mydat{i}), nanvec(i) = true; else nanvec(i) = false; end
end
mydat(nanvec) = {'NaN'};

% try to replace empty with 'NaN'. We will convert to numeric later.
wh_empty = cellfun(@isempty, mydat);
mydat(wh_empty) = {'NaN'};

% sometimes will read in mixed format: char and numeric. 
% if so, we need to convert ALL to char.  
% all-numeric is handled above, so here convert if mixed.
is_numeric = cellfun(@isnumeric, mydat);
for i = 1:length(mydat)
    if is_numeric(i)
        mydat{i} = num2str(mydat{i});
    end
end
% mydat{is_numeric} = cellfun(@num2str, mydat(is_numeric));

% try to convert to numeric
%     numdat = cellfun(@str2num, mydat);
numdat = cellfun(@str2num, mydat, 'UniformOutput', false);

% if any empty after replacing blanks, this is text (leave as cell)
istext = any(cellfun(@isempty, numdat));

if istext
    % Leave as cell array, but make NaN empty char
    mydat(wh_empty) = {' '};

else  % numeric
    mydat = cat(1, numdat{:});

end

end % process_data_column





%
%
%
% function data_table = process_sheet(data, textdata, has_header_row)
% % Create a table object from separate numeric data (data) and text data (textdata).
%
% % % Variable names
% if has_header_row
%
%     varnames = textdata(1, :);
%
% else
%     varnames = cell(1, size(textdata, 2));  % textdata only seems to be full-width
%
%     for i = 1:length(varnames)
%         varnames{i} = sprintf('Var_%03d', i);
%     end
%
% end
%
% varnames = matlab.lang.makeValidName(varnames);
%
% % % Get vector of text vs. numeric
% % --------------------
%
%
% % % All data
% % --------------------
%
% data_struct = struct();
%
% numeric_counter = 1;        % Numeric array has reduced number of columns
%
% for i = 1:length(varnames)
%
%     fname = varnames{i};
%
%     % figure out if is text
%     mytext = textdata((has_header_row+1):end, i); % skip first row if headers.
%     istext = ~all(cellfun(@isempty, mytext));
%
%     if istext
%
%         data_struct.(fname) = mytext;
%
%     else  % numeric
%
%         data_struct.(fname) = data(:, numeric_counter);
%         numeric_counter = numeric_counter + 1;
%
%     end
%
% end
%
% end
%
