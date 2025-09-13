function bidsMapper = extract_bids_info(bidsFilename)
    % extract_bids_info: Extracts key-value pairs from a BIDS-formatted filename
    %
    % This function parses a filename formatted according to the Brain Imaging Data
    % Structure (BIDS) standard and extracts the key-value pairs embedded within it.
    % The function returns a containers.Map object where each key is a string
    % corresponding to a BIDS key, and each value is the associated value from the filename.
    %
    % Inputs:
    %   bidsFilename - A string or character array representing a BIDS-formatted filename.
    %
    % Outputs:
    %   bidsMapper - A containers.Map object containing key-value pairs extracted
    %                from the BIDS filename.
    %
    % Example:
    %   filename = 'sub-01_ses-02_task-rest_bold.nii';
    %   bidsInfo = extract_bids_info(filename);
    %   % bidsInfo is now a map with keys 'sub', 'ses', 'task', and the corresponding values.
    %
    % Note:
    %   The function expects a properly BIDS-formatted filename. Non-conforming filenames
    %   may result in incorrect or incomplete extraction of key-value pairs.
    %
    % See also: containers.Map, strsplit

    % By Michael Sun, Ph.D. 02/09/2024

    % Check if the input is a string or char array
    if ~(ischar(bidsFilename) || isstring(bidsFilename))
        error('Input must be a string or char array');
    end

    % Split the filename by underscores
    keyValuePairs = strsplit(bidsFilename, '_');

    % Initialize an empty containers.Map
    bidsMapper = containers.Map;

    % Loop through each key-value pair
    for i = 1:length(keyValuePairs)
        % Further split each pair by hyphen
        kv = strsplit(keyValuePairs{i}, '-');

        % Check if the pair is valid (contains a key and a value)
        if length(kv) == 2
            key = kv{1};
            value = kv{2};

            % Add the pair to the map
            bidsMapper(key) = value;
        end
    end
end
