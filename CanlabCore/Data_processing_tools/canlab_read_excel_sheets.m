function data = canlab_read_excel_sheets(filename, name_of_example_header_variable)
% Read an Excel .xls with multiple sheets into a series of table objects,
% auto-detecting the number of header rows.
%
% - Will not work for wide tables with >26*2 columns, without adjustment

sheetList = sheetnames(filename);
data = struct();

for i = 1:length(sheetList)

    sheetName = sheetList{i};

    fprintf('Reading %s: sheet %3.0f of %3.0f\n', sheetName, i, length(sheetList))

    % Read the first 10 rows to locate the header row (adjust range if needed)
    raw = readcell(filename, 'Sheet', sheetName); % , 'Range', 'A1:Z10');
    
    % Search for the row that contains the string 'atlas_names'
    headerRow = [];
    for r = 1:size(raw,1)
        if any(strcmp(raw(r,:), name_of_example_header_variable))
            headerRow = r;
            break;
        end
    end

    if isempty(headerRow)
        % Default if header row not found

        data.(fieldName) = T;
        continue                % to next sheet
        %headerRow = 3;
    end
    
    % Extract the header row data from the read cell array
    headerRowData = raw(headerRow, :);

    % Remove missing values
    mfun = @(x) isa(x, 'missing');
    ismissing = cellfun(mfun, headerRowData);
    headerRowData(ismissing) = [];

    % Determine the last non-empty column in the header row
    lastColIndex = find(~cellfun(@(x) isempty(x) || (ischar(x) && all(isspace(x))), headerRowData), 1, 'last');
    if isempty(lastColIndex)
        error('No header data found in sheet %s.', sheetName);
    end
    % Convert the column index to its Excel letter equivalent (assuming columns start at A)
    if lastColIndex <= 25
        lastColLetter = char(double('A') + lastColIndex - 1);
    else
        % more columns, add double letter
        lastColLetter = ['A' char(double('A') + lastColIndex - 27)];
    end

    % Create import options for the sheet and update them
    opts = detectImportOptions(filename, 'Sheet', sheetName);

    % Set the header (variable names) range exactly over the non-empty header cells
    opts.VariableNamesRange = sprintf('A%d:%s%d', headerRow, lastColLetter, headerRow);

    % Tell MATLAB where the data begins (i.e. the row after the header row)
    opts.DataRange = sprintf('A%d', headerRow+1);
    
    % Also, update the VariableNames property to match exactly the non-empty header cells
    opts.VariableNames = headerRowData(1, 1:lastColIndex);
    
    % Read the table using the updated options
    T = readtable(filename, opts);
    
    % Store the table in a structure field with a valid MATLAB name
    fieldName = matlab.lang.makeValidName(sheetName);
    data.(fieldName) = T;

end % loop


end % function