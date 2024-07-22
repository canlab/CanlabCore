function gen_design_files(bids_events_tsv, formats, varargin)
    % gen_design_files generates design files for specified brain imaging analysis software
    % from a BIDS events .tsv file.
    %
    % Usage:
    %   gen_design_files(bids_events_tsv, 'SPM')
    %   gen_design_files(bids_events_tsv, {'FSL', 'SPM'}, 'FSLModulator', 'reaction_time', 'outdir', '/path/to/output')
    %   gen_design_files(bids_events_tsv, {'AFNI', 'FSL'}, 'outdir', '/path/to/output')
    %
    % Inputs:
    %   bids_events_tsv - String, path to the BIDS events .tsv file.
    %   formats - A char or cell array of chars specifying which design file formats to generate.
    %             Acceptable formats are 'SPM', 'FSL', 'AFNI'.
    %   varargin - Variable input arguments:
    %       'outdir' - A string specifying the output directory path.
    %       'FSLModulator' - A string specifying a column in your bids events.tsv file
    %                       to use as a parametric modulator. If not specified, defaults to 1 for all events.
    %
    % Outputs:
    %   The function will write design files to the specified output directory according to the formats specified.
    %
    %  05/07/2024 Michael Sun, Ph.D.
    
    
    % Initialize the parser
    p = inputParser;
    addRequired(p, 'bids_events_tsv', @ischar);
    addRequired(p, 'formats', @(x) ischar(x) || iscellstr(x));
    addParameter(p, 'outdir', pwd, @ischar);  % Default output directory
    addParameter(p, 'FSLModulator', '', @ischar);  % Optional modulator column

    % Parse inputs
    parse(p, bids_events_tsv, formats, varargin{:});

    % Check formats and set flags
    if ischar(formats)
        formats = {formats};  % Convert to cell array if only one format is provided
    end

    generateSPM = ismember('SPM', formats);
    generateFSL = ismember('FSL', formats);
    generateAFNI = ismember('AFNI', formats);

    outdir = p.Results.outdir;
    fslModulator = p.Results.FSLModulator;

    % Check the BIDS events file existence
    if ~exist(bids_events_tsv, 'file')
        error('The specified BIDS events file does not exist.');
    end

    % Load the BIDS events file
    opts = detectImportOptions(bids_events_tsv, 'FileType', 'text');
    bids_data = readtable(bids_events_tsv, opts);

    % Ensure the output directory exists
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % Generate design files based on the formats requested
    if generateSPM
        createSPMDesign(bids_data, outdir);
        fprintf('SPM design files created in %s\n', outdir);
    end

    if generateFSL
        if ~isempty(fslModulator)
            createFSLDesign(bids_data, outdir, fslModulator);
        else
            createFSLDesign(bids_data, outdir);
        end
        fprintf('FSL design files created in %s\n', outdir);
    end

    if generateAFNI
        createAFNIDesign(bids_data, outdir);
        fprintf('AFNI design files created in %s\n', outdir);
    end
end

function createSPMDesign(bids_data, outdir)
    % Extract unique trial types
    unique_types = unique(bids_data.trial_type);
    
    % Create one design file per trial type
    for i = 1:length(unique_types)
        type = unique_types{i};
        
        % Filter rows for the current trial type
        idx = strcmp(bids_data.trial_type, type);
        
        % Collect the names, onsets, and durations
        name = {type};
        onset = {bids_data.onset(idx)};
        duration = {bids_data.duration(idx)};
        
        % Define the output filename
        filename = fullfile(outdir, sprintf('%s.mat', type));
        
        % Save the design file
        save(filename, 'name', 'onset', 'duration');
        fprintf('SPM design file for %s created at %s\n', type, filename);
    end
end

function createFSLDesign(bids_data, outdir, modulatorColumnName)
    % Create FSL design files for each unique trial type, allowing for a specified parametric modulator column.
    %
    % Inputs:
    %   bids_data - Table containing the BIDS event data.
    %   outdir - Directory where the design files will be saved.
    %   modulatorColumnName - (Optional) Name of the column to use for parametric modulators.
    %                         If not provided, all modulators are set to 1.
    
    % Check if a modulator column name is provided and valid
    if nargin < 3 || isempty(modulatorColumnName) || ~ismember(modulatorColumnName, bids_data.Properties.VariableNames)
        modulatorColumnName = ''; % Use default value
        defaultModulatorValue = 1;
    end
    
    % Extract unique trial types
    unique_types = unique(bids_data.trial_type);
    
    % Create one design file per trial type
    for i = 1:length(unique_types)
        type = unique_types{i};
        
        % Filter rows for the current trial type
        idx = strcmp(bids_data.trial_type, type);
        
        % Prepare the data
        onsets = bids_data.onset(idx);
        durations = bids_data.duration(idx);
        if isempty(modulatorColumnName)
            modulators = repmat(defaultModulatorValue, sum(idx), 1); % Default modulators
        else
            modulators = bids_data.(modulatorColumnName)(idx);
        end
        
        % Define the output filename
        filename = fullfile(outdir, sprintf('%s.txt', type));
        
        % Write the data to a file
        data_to_write = [onsets, durations, modulators];
        writematrix(data_to_write, filename, 'Delimiter', ' ');
        fprintf('FSL design file for %s created at %s\n', type, filename);
    end
end

function createAFNIDesign(bids_data, outdir)
    % Create AFNI design files for each unique trial type, formatted for use with AFNI's tools.
    %
    % Inputs:
    %   bids_data - Table containing the BIDS event data.
    %   outdir - Directory where the design files will be saved.

    % Ensure the output directory exists
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % Extract unique trial types
    unique_types = unique(bids_data.trial_type);
    
    % Create one design file per trial type
    for i = 1:length(unique_types)
        type = unique_types{i};
        
        % Filter rows for the current trial type
        idx = strcmp(bids_data.trial_type, type);
        
        % Prepare the data
        onsets = bids_data.onset(idx);
        
        % Convert onsets to a string of space-separated values
        onset_str = strjoin(arrayfun(@num2str, onsets, 'UniformOutput', false), ' ');
        
        % Define the output filename
        filename = fullfile(outdir, sprintf('%s.1D', type));
        
        % Write the data to a file
        fid = fopen(filename, 'wt');
        if fid == -1
            error('Cannot open file %s for writing.', filename);
        end
        fprintf(fid, '%s\n', onset_str);
        fclose(fid);
        
        fprintf('AFNI design file for %s created at %s\n', type, filename);
    end
end