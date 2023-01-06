function canlab_prep_bidsdir(bidsdir, varargin)
    % This function 'canlabulates' a directory, converting a BIDS directory
    % to be canlab_glm_subject_levels compatible such that first-level analyses can be done. 
    % 
    % required argument:
    % 1. bidsdir: full filepath to the BIDS-directory of your study. 
    % e.g., 'F:\Dropbox (Dartmouth College)\Tor Pinel datasets\Pinel_localizer\data\pinel_localizer_Dartmouth_S2019' 
    % 
    % optional varargin:
    % 2. funcstring: the functional filename identifiers (with wildcards) to pull from every subject. 
    % default: '*preproc*.nii.gz' 
    % 3. noisestring: the confound filename identifiers (with wildcards) to pull from every subject. 
    % default: '*confound*.tsv'
    %
    % Function by Michael Sun, Ph.D. 11/29/2022

    % Process arguments:
    if nargin==1
        funcstring='*preproc*.nii.gz';
        noisestring='*confound*.tsv';
    else
        if nargin==2
            funcstring=varargin{1};
        elseif nargin==3
            funcstring=varargin{1};
            noisestring=varargin{2};
        end
    end
    
    % Start:
    if contains(bidsdir, 'func')
        funcs=dir(fullfile(bidsdir, funcstring));   % Identify the functional scans.
        noise=dir(fullfile(bidsdir, noisestring));  % Identify the confound files.
    elseif contains(bidsdir, 'ses')
        funcs=dir(fullfile(bidsdir, 'func', funcstring));   % Identify the functional scans.
        noise=dir(fullfile(bidsdir, 'func', noisestring));  % Identify the confound files.
    else
        funcs=dir(fullfile(bidsdir, 'sub-*', '**', 'func', funcstring));   % Identify the functional scans.
        noise=dir(fullfile(bidsdir, 'sub-*', '**', 'func', noisestring));  % Identify the confound files.
    end

    for f = 1:numel(funcs)
        img=cell2mat(erase(fullfile({funcs(f).folder}, {funcs(f).name}), '.gz'));
        rundir=cell2mat(fullfile(funcs(f).folder, strcat('run-', extractBetween(img, 'run-', '_'))));
        % create run folders if they don't exist yet.
        if ~isdir(rundir)
            mkdir(rundir);
        end
        if ~isfile(img)   % Check if .nii exists
            gunzip(fullfile({funcs(f).folder}, {funcs(f).name})); % Unzip any nii.gz file so that the .nii is revealed; not needed after being run once.
        end
        % Copy the requisite func file to the new folder
        copyfile(img, rundir);
    end

    for n = 1:numel(noise)
        % Copy the requisite noise file to the new folder
        rundir=cell2mat(fullfile(noise(n).folder, strcat('run-', extractBetween(noise(n).name, 'run-', '_'))));
        copyfile(fullfile(noise(n).folder, noise(n).name), rundir);
    end

end