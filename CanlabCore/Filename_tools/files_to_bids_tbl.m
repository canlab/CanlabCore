function t = files_to_bids_tbl(files)
    % filestoBIDStbl: Create a table from BIDS-formatted filenames
    %
    % This function takes an array of BIDS-formatted filenames and converts them
    % into a table. Each row corresponds to a file, and columns are dynamically
    % generated based on the keys present in the filenames.
    %
    % Inputs:
    %   files - A cell array or string array of BIDS-formatted filenames.
    %
    % Outputs:
    %   t - A table containing the extracted BIDS information.
    %
    % Example:
    %   files = ["sub-01_ses-02_task-rest_bold.nii", "sub-02_ses-01_task-run_bold.nii"];
    %   bidsTable = filestoBIDStbl(files);
    %
    % See also: extractBIDSinfo, table
    % Michael Sun, Ph.D. 02/09/2024

    % Gather all unique keys from all files
    allKeys = {};
    for i = 1:numel(files)
        bidsInfo = extract_bids_info(files{i});
        allKeys = union(allKeys, keys(bidsInfo));
    end

    % Initialize the table with filenames and dynamic keys
    variableNames = ['filename', allKeys'];
    t = table('Size', [numel(files), numel(variableNames)], ...
              'VariableTypes', repmat({'string'}, 1, numel(variableNames)), ...
              'VariableNames', variableNames);

    % Process each file and populate the table
    for i = 1:numel(files)
        filename = files{i};
        t.filename(i) = filename;
        bidsInfo = extract_bids_info(filename);

        for key = 1:numel(allKeys)
            if isKey(bidsInfo, allKeys{key})
                t{i, allKeys{key}} = {bidsInfo(allKeys{key})};
            else
                t{i, allKeys{key}} = ""; % Fill with empty string if key is not present
            end
        end
    end
end
