function DSGN = promptDSGN()
    % Build a DSGN structure to be passed into canlab_glm_subject_levels(DSGN) based on a series of user prompts and validate it.
    % Michael Sun, Ph.D.

    % Initialize the DSGN structure
    DSGN = struct();

    % Prompt the user for metadata
    DSGN.metadata = input('Enter metadata (optional): ', 's');

    % Prompt the user for the model directory
    DSGN.modeldir = input('Enter model directory: ', 's');

    % Prompt the user for the subjects' directories
    DSGN.subjects = input('Enter subjects directories (cell array): ');

    % Prompt the user for the functional file names
    DSGN.funcnames = input('Enter functional file names (cell array): ');

    % Prompt the user for optional fields
    DSGN.allowmissingfunc = input('Allow missing functional files? (true/false, optional): ');

    % Check if the field is optional and allow skipping
    if isempty(DSGN.allowmissingfunc)
        DSGN = rmfield(DSGN, 'allowmissingfunc');
    end

    DSGN.concatenation = input('Enter concatenation information (cell array, optional): ');

    if isempty(DSGN.concatenation)
        DSGN = rmfield(DSGN, 'concatenation');
    end

    DSGN.customrunintercepts = input('Enter custom run intercepts (cell array, optional): ');

    if isempty(DSGN.customrunintercepts)
        DSGN = rmfield(DSGN, 'customrunintercepts');
    end

    DSGN.ar1 = input('Use AR(1) to model serial correlations? (true/false, optional): ');

    if isempty(DSGN.ar1)
        DSGN = rmfield(DSGN, 'ar1');
    end

    DSGN.fast = input('Use FAST autoregressive model? (true/false, optional): ');

    if isempty(DSGN.fast)
        DSGN = rmfield(DSGN, 'fast');
    end

    DSGN.notimemod = input('Turn off time modulation of conditions? (true/false, optional): ');

    if isempty(DSGN.notimemod)
        DSGN = rmfield(DSGN, 'notimemod');
    end

    DSGN.singletrials = input('Enter single trial conditions (cell array, optional): ');

    if isempty(DSGN.singletrials)
        DSGN = rmfield(DSGN, 'singletrials');
    end

    DSGN.singletrialsall = input('Set single trials to true for all conditions? (true/false, optional): ');

    if isempty(DSGN.singletrialsall)
        DSGN = rmfield(DSGN, 'singletrialsall');
    end

    DSGN.modelingfilesdir = input('Enter modeling files directory (optional): ', 's');

    if isempty(DSGN.modelingfilesdir)
        DSGN = rmfield(DSGN, 'modelingfilesdir');
    end

    DSGN.allowemptycond = input('Allow empty conditions? (true/false, optional): ');

    if isempty(DSGN.allowemptycond)
        DSGN = rmfield(DSGN, 'allowemptycond');
    end

    DSGN.allowmissingcondfiles = input('Allow missing condition files? (true/false, optional): ');

    if isempty(DSGN.allowmissingcondfiles)
        DSGN = rmfield(DSGN, 'allowmissingcondfiles');
    end

    % Prompt the user for additional optional fields with examples
    DSGN.multireg = input('Enter multiple regressors file (optional, e.g., ''noise_model1''): ', 's');

    if isempty(DSGN.multireg)
        DSGN = rmfield(DSGN, 'multireg');
    end

    DSGN.multiregbehav = input('Enter multiple regressors file (behavioral) (optional): ', 's');

    if isempty(DSGN.multiregbehav)
        DSGN = rmfield(DSGN, 'multiregbehav');
    end

    DSGN.contrasts = input('Enter contrasts (cell array, optional): ');

    if isempty(DSGN.contrasts)
        DSGN = rmfield(DSGN, 'contrasts');
    end

    DSGN.contrastnames = input('Enter contrast names (cell array, optional): ');

    if isempty(DSGN.contrastnames)
        DSGN = rmfield(DSGN, 'contrastnames');
    end

    DSGN.contrastweights = input('Enter contrast weights (cell array, optional): ');

    if isempty(DSGN.contrastweights)
        DSGN = rmfield(DSGN, 'contrastweights');
    end

    DSGN.regmatching = input('Enter regular expression matching (optional): ', 's');

    if isempty(DSGN.regmatching)
        DSGN = rmfield(DSGN, 'regmatching');
    end

    DSGN.defaultsuffix = input('Enter default suffix (optional): ', 's');

    if isempty(DSGN.defaultsuffix)
        DSGN = rmfield(DSGN, 'defaultsuffix');
    end

    DSGN.noscale = input('Use noscale option? (true/false, optional): ');

    if isempty(DSGN.noscale)
        DSGN = rmfield(DSGN, 'noscale');
    end

    % Validate the provided inputs
    validateDSGN(DSGN);

end

function validateDSGN(DSGN)

    % Validation for DSGN.metadata
    if ~isfield(DSGN, 'metadata') || ~ischar(DSGN.metadata)
        error('Invalid or missing value for DSGN.metadata.');
    end
    
    % Validation for DSGN.modeldir
    if ~isfield(DSGN, 'modeldir') || ~ischar(DSGN.modeldir)
        error('Invalid or missing value for DSGN.modeldir.');
    end

    % Check if model directory exists
    if ~isfolder(DSGN.modeldir)
        error('Invalid model directory');
    end
    
    % Validation for DSGN.subjects
    if ~isfield(DSGN, 'subjects') || ~iscell(DSGN.subjects) || isempty(DSGN.subjects)
        error('Invalid or missing value for DSGN.subjects.');
    end

    % Check if subjects directories are valid
    for i = 1:numel(DSGN.subjects)
        if ~isfolder(DSGN.subjects{i})
            error('Invalid subjects directory at index %d', i);
        end
    end

    
    % Validation for DSGN.funcnames
    if ~isfield(DSGN, 'funcnames') || ~iscell(DSGN.funcnames) || isempty(DSGN.funcnames)
        error('Invalid or missing value for DSGN.funcnames.');
    end

    % Check if functional files exist
    for i = 1:numel(DSGN.funcnames)
        filenames = dir(fullfile(DSGN.subjects{1}, DSGN.funcnames{i}));
        if isempty(filenames)
            error('Functional file not found for session %d', i);
        end
    end


    
    % Validation for DSGN.tr
    if ~isfield(DSGN, 'tr') || ~isnumeric(DSGN.tr) || numel(DSGN.tr) ~= 1 || isnan(DSGN.tr) || DSGN.tr <= 0
        error('Invalid or missing value for DSGN.tr.');
    end
    
    % Validation for DSGN.hpf
    if ~isfield(DSGN, 'hpf') || ~isnumeric(DSGN.hpf) || numel(DSGN.hpf) ~= 1 || isnan(DSGN.hpf) || DSGN.hpf < 0
        error('Invalid or missing value for DSGN.hpf.');
    end
    
    % Validation for DSGN.fmri_t
    if ~isfield(DSGN, 'fmri_t') || ~isnumeric(DSGN.fmri_t) || numel(DSGN.fmri_t) ~= 1 || isnan(DSGN.fmri_t) || DSGN.fmri_t <= 0
        error('Invalid or missing value for DSGN.fmri_t.');
    end
    
    % Validation for DSGN.fmri_t0
    if ~isfield(DSGN, 'fmri_t0') || ~isnumeric(DSGN.fmri_t0) || numel(DSGN.fmri_t0) ~= 1 || isnan(DSGN.fmri_t0) || DSGN.fmri_t0 <= 0
        error('Invalid or missing value for DSGN.fmri_t0.');
    end
    
    % Validation for DSGN.conditions
    if ~isfield(DSGN, 'conditions') || ~iscell(DSGN.conditions) || isempty(DSGN.conditions)
        error('Invalid or missing value for DSGN.conditions.');
    end
    
    % Validation for DSGN.allowmissingfunc (optional)
    if isfield(DSGN, 'allowmissingfunc') && (~islogical(DSGN.allowmissingfunc) || numel(DSGN.allowmissingfunc) ~= 1)
        error('Invalid value for DSGN.allowmissingfunc.');
    end
    
    % Validation for DSGN.concatenation (optional)
    if isfield(DSGN, 'concatenation') && (~iscell(DSGN.concatenation) || isempty(DSGN.concatenation))
        error('Invalid value for DSGN.concatenation.');
        % Additional validation logic specific to DSGN.concatenation can be added here
    end
    
    % Validation for DSGN.customrunintercepts (optional)
    if isfield(DSGN, 'customrunintercepts') && (~iscell(DSGN.customrunintercepts) || isempty(DSGN.customrunintercepts))
        error('Invalid value for DSGN.customrunintercepts.');
        % Additional validation logic specific to DSGN.customrunintercepts can be added here
    end
    
    % Validation for DSGN.ar1 (optional)
    if isfield(DSGN, 'ar1') && (~islogical(DSGN.ar1) || numel(DSGN.ar1) ~= 1)
        error('Invalid value for DSGN.ar1.');
    end
    
    % Validation for DSGN.fast (optional)
    if isfield(DSGN, 'fast') && (~islogical(DSGN.fast) || numel(DSGN.fast) ~= 1)
        error('Invalid value for DSGN.fast.');
    end

    % Check if AR(1) and FAST options are mutually exclusive
    if isfield(DSGN, 'ar1') && isfield(DSGN, 'fast')
        if DSGN.ar1 && DSGN.fast
            error('AR(1) and FAST options are mutually exclusive');
        end
    end

    
    % Validation for DSGN.notimemod (optional)
    if isfield(DSGN, 'notimemod') && (~islogical(DSGN.notimemod) || numel(DSGN.notimemod) ~= 1)
        error('Invalid value for DSGN.notimemod.');
    end
    
    % Validation for DSGN.singletrials (optional)
    if isfield(DSGN, 'singletrials') && (~iscell(DSGN.singletrials) || isempty(DSGN.singletrials))
        error('Invalid value for DSGN.singletrials.');
        % Additional validation logic specific to DSGN.singletrials can be added here
    end

    % Validation for DSGN.singletrialsall (optional)
    if isfield(DSGN, 'singletrialsall') && (~islogical(DSGN.singletrialsall) || numel(DSGN.singletrialsall) ~= 1)
        error('Invalid value for DSGN.singletrialsall.');
    end

    % Check if singletrials and singletrialsall options are mutually exclusive
    if isfield(DSGN, 'singletrials') && isfield(DSGN, 'singletrialsall')
        if ~isempty(DSGN.singletrials) && DSGN.singletrialsall
            error('singletrials and singletrialsall options are mutually exclusive');
        end
    end
    
    % Validation for DSGN.modelingfilesdir (optional)
    if isfield(DSGN, 'modelingfilesdir') && (~ischar(DSGN.modelingfilesdir) || isempty(DSGN.modelingfilesdir))
        error('Invalid value for DSGN.modelingfilesdir.');
    end
    
    % Validation for DSGN.allowemptycond (optional)
    if isfield(DSGN, 'allowemptycond') && (~islogical(DSGN.allowemptycond) || numel(DSGN.allowemptycond) ~= 1)
        error('Invalid value for DSGN.allowemptycond.');
    end
    
    % Validation for DSGN.allowmissingcondfiles (optional)
    if isfield(DSGN, 'allowmissingcondfiles') && (~islogical(DSGN.allowmissingcondfiles) || numel(DSGN.allowmissingcondfiles) ~= 1)
        error('Invalid value for DSGN.allowmissingcondfiles.');
    end
    
    % Validation for DSGN.multireg (optional)
    if isfield(DSGN, 'multireg') && (~ischar(DSGN.multireg) || isempty(DSGN.multireg))
        error('Invalid value for DSGN.multireg.');
    end
    
    % Validation for DSGN.multiregbehav (optional)
    if isfield(DSGN, 'multiregbehav') && (~ischar(DSGN.multiregbehav) || isempty(DSGN.multiregbehav))
        error('Invalid value for DSGN.multiregbehav.');
    end
    
    % Validation for DSGN.contrasts (optional)
    if isfield(DSGN, 'contrasts') && (~iscell(DSGN.contrasts) || isempty(DSGN.contrasts))
        error('Invalid value for DSGN.contrasts.');
    end

    % Check if defaultsuffix is provided when regmatching is set to 'regexp'
    if isfield(DSGN, 'regmatching') && strcmp(DSGN.regmatching, 'regexp')
        if ~isfield(DSGN, 'defaultsuffix')
            error('defaultsuffix must be provided when regmatching is set to ''regexp''');
        end
    end


    % Display a success message if all validations pass
    disp('DSGN structure created successfully.');

end