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
    %   canlab_prep_bidsdir(bidsdir, 'func', funcstring, 'noise', noisestring, 'outdir', outdir)
    %
    % Required Argument:
    %   bidsdir       - Full filepath to the BIDS-directory of your study.
    %                   e.g., 'F:\Dropbox (Dartmouth College)\Tor Pinel datasets\Pinel_localizer\data\pinel_localizer_Dartmouth_S2019'
    %
    % Optional Key-Value Arguments:
    %   'func'        - Functional filename identifiers (with wildcards) to pull from 
    %                   every subject. Default: '*preproc*.nii.gz'.
    %
    %   'noise'       - Confound filename identifiers (with wildcards) to pull from 
    %                   every subject. Default: '*confound*.tsv'.
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


    % Additional argument to handle custom output directory

    % Create an input parser object
    p = inputParser;

    % Define default values
    defaultFuncString = '*preproc*.nii.gz';
    defaultNoiseString = '*confound*.tsv';
    defaultOutDir = bidsdir;

    % Add required and optional inputs
    addRequired(p, 'bidsdir', @ischar);
    addParameter(p, 'func', defaultFuncString, @ischar);
    addParameter(p, 'noise', defaultNoiseString, @ischar);
    addParameter(p, 'outdir', defaultOutDir, @ischar);

    % Parse the inputs
    parse(p, bidsdir, varargin{:});

    % Extract the parsed values
    funcstring = p.Results.func;
    noisestring = p.Results.noise;
    output_dir = p.Results.outdir;


    [~, lastdir] = fileparts(bidsdir);    % get the last part of the directory

    % Start:
    if contains(lastdir, 'func')
        funcs=dir(fullfile(bidsdir, funcstring));   % Identify the functional scans.
        noise=dir(fullfile(bidsdir, noisestring));  % Identify the confound files.
    elseif contains(lastdir, 'ses')
        funcs=dir(fullfile(bidsdir, 'func', funcstring));   % Identify the functional scans.
        noise=dir(fullfile(bidsdir, 'func', noisestring));  % Identify the confound files.
    else
        funcs=dir(fullfile(bidsdir, 'sub-*', '**', 'func', funcstring));   % Identify the functional scans.
        noise=dir(fullfile(bidsdir, 'sub-*', '**', 'func', noisestring));  % Identify the confound files.
    end

    h = waitbar(0, 'Processing func files...');
    for f = 1:numel(funcs)
        img = cell2mat(fullfile({funcs(f).folder}, {funcs(f).name}));

        % Recreate the original directory structure in output_dir
        relative_path = fileparts(strrep(img, bidsdir, ''));  % Extract the relative path
        disp(relative_path)
        new_path = fullfile(output_dir, relative_path);       % Combine with output_dir

        disp(new_path)
        if ~isdir(new_path)
            mkdir(new_path);
        end

        rundir = cell2mat(fullfile(new_path, strcat('run-', extractBetween(img, 'run-', '_'))));
        % create run folders if they don't exist yet.
        if ~isdir(rundir)
            mkdir(rundir);
        end

        % Create a symlink for the requisite func file to the new folder
        if ispc  % Check if the system is Windows
            cmd_str = ['cmd.exe /C mklink "' fullfile(rundir, funcs(f).name) '" "' img '"'];
            system(cmd_str);
        else
            system(['ln -s ' img ' ' rundir]);
        end
        % Update waitbar
        waitbar(f / numel(funcs), h, sprintf('Processing func file %d of %d...(%d%%)', f, numel(funcs), round((f/numel(funcs))*100)));
    end
    close(h);

    h = waitbar(0, 'Processing noise files...');
    disp('Func files all symlinked.')

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