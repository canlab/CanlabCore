function canlab_prep_bidsdir(bidsdir, varargin)
    % CANLAB_PREP_BIDSDIR Convert a BIDS directory for canlab_glm_subject_levels compatibility
    %
    % This function 'canlabulates' a directory, making it compatible for 
    % canlab_glm_subject_levels first-level analyses. It creates symbolic links
    % instead of copying files, ensuring no redundant space usage while keeping
    % the original BIDS directory intact. Optionally, you can create a BIDS
    % structure outside the original BIDS directory by specifying 'outdir'.
    %
    % Usage:
    %   canlab_prep_bidsdir(bidsdir, 'datatype', datatypestring, 'noise', noisestring, 'outdir', outdir)
    %
    % Required Argument:
    %   bidsdir       - Full filepath to the BIDS-directory of your study.
    %                   e.g., 'F:\Dropbox (Dartmouth College)\Tor Pinel datasets\Pinel_localizer\data\pinel_localizer_Dartmouth_S2019'
    %
    % Optional Key-Value Arguments:
    %   'datatype'    - Data type ('anat' or 'func') filename identifiers (with wildcards) to pull from 
    %                   every subject. Default: 'func'.
    %   'filename'    - The filename identifiers (with wildcards) based on datatype.
    %                   Default: '*space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz' for 'func' and *T1w.nii.gz' for 'anat'.
    %
    %   'noise'       - Confound filename identifiers (with wildcards) to pull from 
    %                   every subject. Default: '*desc-confounds_timeseries.tsv'.
    %
    %   'outdir'      - Full filepath of the directory where the new structure with symbolic links will be created.
    %                   By default, it's the same as the input BIDS directory.
    %
    % Notes:
    %   - For symbolic link creation on Windows, you must run Matlab as administrator.
    %   - Always ensure you have backup copies of your data before performing directory operations.
    %
    % Function by Michael Sun, Ph.D. 11/29/2022
    % Updated: 10/19/2023

    % Create an input parser object
    p = inputParser;

    % Define default values
    defaultDataType = 'func';
    defaultFuncString = '*space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz';
    defaultAnatString = '*T1w.nii.gz'; % Choose a suitable default for anat
    defaultNoiseString = '*desc-confounds_timeseries.tsv';
    defaultOutDir = bidsdir;

    % Add required and optional inputs
    addRequired(p, 'bidsdir', @ischar);
    addParameter(p, 'datatype', defaultDataType, @(x) any(validatestring(x, {'func', 'anat'})));
    addParameter(p, 'filename', defaultFuncString, @ischar);
    addParameter(p, 'noise', defaultNoiseString, @ischar);
    addParameter(p, 'outdir', defaultOutDir, @ischar);

    % Parse the inputs
    parse(p, bidsdir, varargin{:});

    % Extract the parsed values
    datatype = p.Results.datatype;
    if strcmp(datatype, 'func')
        datastring = p.Results.filename;
    else % anat
        datastring = p.Results.filename;
    end
    noisestring = p.Results.noise;
    output_dir = p.Results.outdir;

    % Override the filename default based on datatype
    if strcmp(datatype, 'anat')
        datastring = defaultAnatString;
    elseif strcmp(datatype, 'func')
        datastring = defaultFuncString;
    else
        error('Unsupported datatype. Must be "func" or "anat".');
    end
    
    if strcmp(p.Results.filename, defaultFuncString) || strcmp(p.Results.filename, defaultAnatString)
        % If the filename is default, then override with datastring
        filename = datastring;
    else
        filename = p.Results.filename;
    end

    [~, lastdir] = fileparts(bidsdir);

    % Start:
    if contains(lastdir, datatype)
        datafiles = dir(fullfile(bidsdir, filename));   % Identify the data files.
        if strcmp(datatype, 'func')
            noise = dir(fullfile(bidsdir, noisestring));      % Identify the confound files.
        end
    elseif contains(lastdir, 'ses')
        datafiles = dir(fullfile(bidsdir, datatype, filename));   % Identify the data files.
        if strcmp(datatype, 'func')
            noise = dir(fullfile(bidsdir, datatype, noisestring));      % Identify the confound files.
        end
    else
        datafiles = dir(fullfile(bidsdir, 'sub-*', '**', datatype, filename));   % Identify the data files.
        if strcmp(datatype, 'func')
            noise = dir(fullfile(bidsdir, 'sub-*', '**', datatype, noisestring));      % Identify the confound files.
        end
    end

    h = waitbar(0, ['Processing ' datatype ' files...']);
    for f = 1:numel(datafiles)
        img = cell2mat(fullfile({datafiles(f).folder}, {datafiles(f).name}));

        % Recreate the original directory structure in output_dir
        relative_path = fileparts(strrep(img, bidsdir, ''));  % Extract the relative path
        disp(relative_path)
        new_path = fullfile(output_dir, relative_path);       % Combine with output_dir

        disp(new_path)
        if ~isdir(new_path)
            mkdir(new_path);
        end

        if strcmp(datatype, 'func')
            rundir = cell2mat(fullfile(new_path, strcat('run-', extractBetween(img, 'run-', '_'))));
            % create run folders if they don't exist yet.
            if ~isdir(rundir)
                mkdir(rundir);
            end
    
            % Create a symlink for the requisite func file to the new folder
            if ispc  % Check if the system is Windows
                cmd_str = ['cmd.exe /C mklink "' fullfile(rundir, datafiles(f).name) '" "' img '"'];
                system(cmd_str);
            else
                system(['ln -s ' img ' ' rundir]);
            end
        else
            % Create a symlink for the requisite func file to the new folder
            if ispc  % Check if the system is Windows
                cmd_str = ['cmd.exe /C mklink "' fullfile(new_path, datafiles(f).name) '" "' img '"'];
                system(cmd_str);
            else
                system(['ln -s ' img ' ' new_path]);
            end

        end
        % Update waitbar
        waitbar(f / numel(datafiles), h, sprintf('Processing func file %d of %d...(%d%%)', f, numel(datafiles), round((f/numel(datafiles))*100)));
    end
    close(h);

    if strcmp(datatype, 'func')
    
        h = waitbar(0, 'Processing noise files...');
        disp(['Func ' datatype ' all symlinked.'])
    
        for n = 1:numel(noise)
            % Recreate the original directory structure in output_dir
            relative_path = fileparts(strrep(img, bidsdir, ''));  % Extract the relative path
            disp(relative_path)
            new_path = fullfile(output_dir, relative_path);       % Combine with output_dir
    
            disp(new_path)
            if ~isdir(new_path)
                mkdir(new_path);
            end
            
            % Create a symlink for the requisite noise file to the new folder
            rundir = cell2mat(fullfile(output_dir, strcat('run-', extractBetween(noise(n).name, 'run-', '_'))));
            if ispc  % Check if the system is Windows
                cmd_str = ['cmd.exe /C mklink "' fullfile(rundir, noise(n).name) '" "' fullfile(noise(n).folder, noise(n).name) '"'];
                system(cmd_str);
            else
                system(['ln -s ' fullfile(noise(n).folder, noise(n).name) ' ' rundir]);
            end
    
            % Update waitbar
            waitbar(n / numel(noise), h, sprintf('Processing noise file %d of %d...(%d%%)', n, numel(noise), round((n/numel(noise))*100)));
    
        end
        disp('Noise files all symlinked. Done.')
    
        close(h);
    end
end