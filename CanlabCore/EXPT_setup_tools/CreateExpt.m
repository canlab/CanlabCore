function EXPT = CreateExpt(Action,varargin)
% EXPT = CreateExpt(Action,varargin)
%
% Creates an EXPT experiment information structure for particular analysis programs
% Checks for needed fields for the particular analysis type and attempts to
% create them if they're not there.
%
% Uses global EXPT structure.
% You must save EXPT.mat yourself after you create it.
%
%
% Tor Wager
% Last modified: Oct 2006; July 2007
%
% Uses:
%
% Used in BrainTools GUI and with other toolboxes
%
% Methods (Actions) are:
% 'fir','hewma', 'parammod' 'robfit' 'extractdx' 'CreateField'
%
% General methods for setting up analysis-specific information structures
% -------------------------------------------------------------------------
% EXPT = CreateExpt('Init');    % check for existing EXPT, load existing or use in-memory
% EXPT = CreateExpt('fir');     % set up FIR whole-brain model info
% EXPT = CreateExpt('hewma');   % set up HEWMA whole-brain model info
% EXPT = CreateExpt('parammod');    % set up parametric modulator model info
% EXPT = CreateExpt('robfit');  % set up Robust Regression analysis info
% EXPT = CreateExpt('extractdx');   % set up FIR estimate ROI data extraction info
%
% Methods for setting up specific fields in the EXPT structure
% -------------------------------------------------------------------------
% EXPT = CreateExpt('CreateField',fieldname,substructure);
% Add a field if it doesn't already exist, in substructure specified by
% string
% Examples are:
% EXPT = CreateExpt('CreateField','subjects');    % create EXPT.subjects
% EXPT = CreateExpt('CreateField','P','SNPM');  % create EXPT.SNPM.P
% EXPT = CreateExpt('CreateField','dxbetas','FIR');  % get hrf timecourse
%                                                   image names
%
% EXPT = CreateExpt('extractdx');  % extraction of FIR timecourse setup
% EXPT = CreateExpt('CreateField','im_files');
%
% tor wager, last update Nov 06


try
    EXPT = evalin('base', 'EXPT');
catch
    EXPT = [];
end



    if exist('EXPT','var') && isempty(EXPT) 
        disp('EXPT is empty.  Assuming this is new EXPT (if you wanted to use a saved EXPT, you may have to load it.).');
    end

    %-----------------------------functions-called------------------------
    %
    %-----------------------------functions-called------------------------


    %-Format arguments
    %-----------------------------------------------------------------------
    if nargin == 0, Action='Init'; end

    switch lower(Action)


        %======================================================================
        %
        %
        % Methods of general utility
        % These are called directly from the command line
        % Or are called by the 'batch' analysis methods below
        %
        %======================================================================

        case lower('Init')
            % ---------------------------------------------------------------------

            % load EXPT
            if isempty(EXPT) && ~(exist('EXPT.mat', 'file'))
                disp('You need to create a variable that has information about the experiment (EXPT).');
                disp('This variable is called EXPT, and is stored in the main study directory in EXPT.mat')
                disp('No EXPT file in current directory. ');
                fprintf(1,'\n')
                
            elseif ~isempty(EXPT)
                disp('Using EXPT already in memory.');
                
            elseif exist('EXPT.mat', 'file')
                disp('loading EXPT.mat from file in current directory'); load EXPT;

            end


        case lower('CreateField')
            % ---------------------------------------------------------------------
            % creates a field with a specific name
            switch length(varargin)
                case 0, error('Must specify field name.');
                case 1, EXPT = CreateField(EXPT, varargin{1});
                case 2, EXPT = CreateField(EXPT, varargin{1},varargin{2});
            end



            %======================================================================
            %
            %
            % Analysis-specific methods
            % These are 'batch' setups for specific analysis toolboxes
            % in the SCANlab tool set.
            %
            %======================================================================

            % ---------------------------------------------------------------------
        case {'fir','hewma'}
            % ---------------------------------------------------------------------

            show_banner(Action);
            EXPT = CreateField(EXPT, 'subjects');              % List of subject directories
            EXPT = CreateField(EXPT, 'im_files','FILES');    % Indiv. subject image names
            EXPT = CreateField(EXPT, 'TR');                  % Scanning TR
            EXPT = CreateField(EXPT, 'HP','FIR');            % High-pass filter
            EXPT = CreateField(EXPT, 'nruns','FIR');
            % nruns = number of sessions (intercepts) OR num. of images in each sess,
            %           e.g., [169 169 172]

            EXPT = CreateField(EXPT, 'numframes','FIR');
            % numframes = number of beta images per event type, e.g., [20 20 20] for FIR
            %           model

            EXPT = CreateField(EXPT, 'smoothlen','FIR');     % FIR beta-series smoothing (time to zero in imgs)

            EXPT = CreateField(EXPT, 'mask');



            % ---------------------------------------------------------------------
        case lower('parammod')
            % ---------------------------------------------------------------------
            disp('In progress.')
            show_banner('ParamMod');
            CreateExpt('FIR')



            % ---------------------------------------------------------------------
        case lower('robfit')
            % ---------------------------------------------------------------------
            show_banner('robfit');
            if isfield(EXPT, 'subjects')
                cnew = input('Found EXPT.subjects.  Re-create? (1/0) ');
                if cnew, EXPT = rmfield(EXPT,'subjects'); end
            end
            EXPT = CreateField(EXPT, 'subjects');               % List of subject directories
            if isfield(EXPT,'SNPM') && isfield(EXPT.SNPM,'P'), EXPT.SNPM = rmfield(EXPT.SNPM,'P'); end
            if isfield(EXPT,'SNPM') && isfield(EXPT.SNPM,'connames'),EXPT.SNPM = rmfield(EXPT.SNPM,'connames'); end
            if isfield(EXPT,'mask'),EXPT = rmfield(EXPT,'mask'); end

            EXPT = CreateField(EXPT, 'P','SNPM');               % images
            EXPT = CreateField(EXPT, 'connames','SNPM');        % Contrast names
            EXPT = CreateField(EXPT, 'mask');



            % ---------------------------------------------------------------------
        case lower('extractdx')
            % ---------------------------------------------------------------------
            % set up required information to extract FIR model H/T/W timecourses
            % into chosen regions of interest (and plot).
            show_banner('extractdx');

            if ~isfield(EXPT.FIR,'dxnames'), error('Please enter cell array of FIR event names in EXPT.FIR.dxnames'); end


            if ~isfield(EXPT.FIR,'TR'),
                if isfield(EXPT,'TR')
                    EXPT.FIR.TR = EXPT.TR;
                else
                    EXPT.FIR.TR = input('Please enter TR in s: ');
                end
            end

            if ~isfield(EXPT.FIR,'baseline'), EXPT.FIR.baseline = input('Please enter indices of baseline timepoints or []: '); end

            if ~isfield(EXPT.FIR,'regsofinterest')
                disp(EXPT.FIR.dxnames)
                wh = input('Enter vector of conditions to plot: ');
                EXPT.FIR.regsofinterest = wh;
            end

            if ~isfield(EXPT.FIR,'mcol')
                mcol = input('Enter vector of colors, e.g. {''ro-'' ''gd'' etc.}: ');
                EXPT.FIR.mcol = mcol;
            end

            if ~isfield(EXPT.FIR,'indiv')
                indiv = input('Enter 1 for breakdown by indiv diffs or 0 for not.');
                EXPT.FIR.indiv = indiv;
                if indiv
                    if ~isfield(EXPT,'beh'), error('NO EXPT.beh vector for indiv diffs -- please enter.'); end
                end
            else
                % do nothing; already exists
                % indiv = EXPT.FIR.indiv;
            end

            EXPT = CreateExpt('CreateField','dxbetas','FIR');




        otherwise
            % ---------------------------------------------------------------------
            error('Unknown action string')

            % ---------------------------------------------------------------------



    end

    % return EXPT to base workspace
    assignin('base', 'EXPT', EXPT);
    
    return









    %=======================================================================
    %=======================================================================

    % sub-functions


    %=======================================================================
    %=======================================================================



function show_banner(str)

    % banner
    fprintf(1,'\n* ======================================================================== *')
    fprintf(1,'\n*                                                                          *')
    fprintf(1,'\n*                 CreateExpt: batch data structure setup                   *')
    fprintf(1,'\n*                                                                          *')
    fprintf(1,'\n*                        Tor Wager, version July 2006                      *')
    fprintf(1,'\n*                                                                          *')
    fprintf(1,'\n* This function creates a data structure called EXPT that stores info      *')
    fprintf(1,'\n* about your experimental designs and parameters.  This structure can be   *')
    fprintf(1,'\n* saved in a file called EXPT.mat and may be useful for reference.         *')
    fprintf(1,'\n* Many analysis toolboxes refer to various fields of EXPT to get filenames *')
    fprintf(1,'\n* and other info needed for analysis.                                      *')
    fprintf(1,'\n* CreateExpt is a utility tool to help you set this up, but you can create *')
    fprintf(1,'\n* or edit any of the fields manually.                                      *')
    fprintf(1,'\n*                                                                          *')
    fprintf(1,'\n* This process is sort of the ''weakest link''   in the analyses at present  *')
    fprintf(1,'\n* Please contact me if you run into bugs.                                  *')
    fprintf(1,'\n*                                                                          *')
    fprintf(1,'\n* ======================================================================== *\n')
    fprintf(1,'Setting up analysis for: %s\n',str);
    fprintf(1,'Required fields that CreateExpt will attempt to create are: \n');

    switch lower(str)

        case {'fir','hewma'}
            reqfields = {'EXPT.subjects' 'EXPT.FILES.im_files' 'EXPT.TR' 'HP' 'EXPT.FIR.nruns' 'EXPT.FIR.numframes' 'EXPT.FIR.smoothlen' 'EXPT.mask'};
            explan = {'List of subject directory names' ...
                'Timeseries image files (full preprocessed data) for each subject in cell array' ...
                'Repetition time (TR) of expt in sec' ...
                'High pass filter cutoff in sec' ...
                'Vector of number of images in each run, e.g., [160 160 180] for a 3-run expt.' ...
                'Number of FIR beta images per event type' ...
                'Exponential smoothing kernel length for FIR estimates (0 or integer)' ...
                'Name of mask image containing in-analysis voxels.'};

        case {'robfit'}
            reqfields = {'EXPT.subjects' 'EXPT.SNPM.P' 'EXPT.SNPM.connames' 'EXPT.SNPM.connums' 'EXPT.mask'};
            explan = {'List of subject directory names' ...
                'Cell array; each cell is a string matrix of image file names, one per subject; one cell per contrast' ...
                'String matrix of names of each contrast, corresponding to cells in EXPT.SNPM.P' ...
                'Numbers of contrast images,corresponding to cells in EXPT.SNPM.P; determines output directory names' ... 
                'Name of mask image containing in-analysis voxels.'};


        case {'extractdx'}

            reqfields = {'EXPT.FIR.dxnames' 'EXPT.FIR.TR' 'EXPT.FIR.baseline' 'EXPT.FIR.regsofinterest' 'EXPT.FIR.mcol' 'EXPT.FIR.indiv' 'EXPT.FIR.dxbetas'};
            explan = {'Cell array of names of conditions for each FIR HRF estimated (create before running!)' ...
                'Repetition time (TR) of expt in sec' ...
                'Baseline timepoints each each FIR HRF (subtracted from HRF) ' ...
                'Indices of which conditions in EXPT.FIR.dxnames are of interest to plot' ...
                'Cell array of colors for each condition of interest' ...
                'Flag to separate plots by individual difference variable (1 or 0)' ...
                'Cell array of str matrices of all images dx* that contain FIR beta estimates, one cell per condition, Ss are rows' };




        otherwise, warning(['Unknown setup command option:' str]);
    end

    for i = 1:length(reqfields), fprintf(1,'%s\t\t%s\n',reqfields{i},explan{i}); end
    fprintf(1,'\n* ======================================================================== *\n')
    fprintf(1,'\n');

    return




function EXPT = CreateField(EXPT, fname,varargin)
    % CreateField(EXPT, fname,[subfield string])
    %
    % Multi-function function to create fields in EXPT
    %

    %=======================================================================
    % check to see if it exists; if so, leave alone
    %=======================================================================
    str = 'EXPT';
    
    if ~isempty(varargin)
        str = [str '.' varargin{1}];

        isparentfield = isfield(EXPT,varargin{1});
        if ~isparentfield, EXPT.(varargin{1}) = []; end

        targetfield = [varargin{1} '.' fname];  % for image name utility
    else
        targetfield = fname;
    end

    try
        is = isfield(eval(str),fname) && ~isempty(eval([str '.' fname]));
    catch
        error('Trying to create subfield when parent doesn''t exist?  Create parent field manually.');
    end

    if is
        fprintf(1,'%s exists.\n',[str '.' fname]);
        return
    else
        fprintf(1,'Creating field: %s.\n',[str '.' fname]);
    end

    %=======================================================================
    % field definition methods
    %=======================================================================
    switch lower(fname)


        % ---------------------------------------------
        % Names of each subject's directory
        % ---------------------------------------------
        case lower('subjects')
            EXPT = get_expt_subdir(EXPT);

            % ---------------------------------------------
            % Image files for single-subject models (processed imgs)
            % ---------------------------------------------
        case lower('im_files')
            fprintf(1,'Getting individual functional images for analysis from subject directories.\n');
            wcard = input('Enter image wildcard (e.g., sn*img): ','s');
            sdir = input('Enter subdirectories within subject dirs to look in (e.g., scan*, or return for none): ','s');
            EXPT = getfunctnames2(EXPT,wcard,targetfield,sdir);
            fprintf(1,'\n');

            % ---------------------------------------------
            % Parameters needed to run single-subject models
            % ---------------------------------------------
        case lower('TR')
            if ~isempty(varargin)
                EXPT.(varargin{1}).(fname) = input('Enter TR in s: ');
            else
                EXPT.(fname) = input('Enter TR in s: ');
            end

        case lower('HP')
            if ~isempty(varargin)
                EXPT.(varargin{1}).(fname) = input('Enter HP filter in s: ');
            else
                EXPT.(fname) = input('Enter HP filter in s: ');
            end

        case lower('nruns')
            quest = sprintf('Enter number of sessions (intercepts) OR num. of images in each sess\ne.g., [169 169 172] for three runs with 169 imgs in first run. : ');
            if ~isempty(varargin)
                EXPT.(varargin{1}).(fname) = input(quest);
            else
                EXPT.(fname) = input(quest);
            end

        case lower('numframes')
            quest = sprintf('Enter number of beta images per event type in FIR, e.g., [20 20 20] : ');
            if ~isempty(varargin)
                EXPT.(varargin{1}).(fname) = input(quest);
            else
                EXPT.(fname) = input(quest);
            end

        case lower('smoothlen')
            quest = sprintf('Enter smoothing length time-to-zero influence (s) for beta FIR estimates (e.g., 6),  : ');
            if ~isempty(varargin)
                EXPT.(varargin{1}).(fname) = input(quest);
            else
                EXPT.(fname) = input(quest);
            end

            % ---------------------------------------------
            % Get an analysis mask name and reslice if necessary
            % ---------------------------------------------
        case lower('mask')
            P = scan_get_files(Inf,'*.img','Select mask image for analysis, or no image for default.',pwd);
            if isempty(P), P = which('scalped_avg152T1_graymatter_smoothed.img'); end
            if ~isempty(varargin)
                EXPT.(varargin{1}).(fname) = P;
            else
                EXPT.(fname) = P;
            end

            try
                % check mask against first image here.
                V1 = spm_vol(EXPT.mask); V1 = V1.mat; V2 = spm_vol(deblank(EXPT.FILES.im_files{1}(1,:))); V2 = V2.mat;
                go = V1 - V2; go = any(go(:));
                if go
                    [P,EXPT.(fname)] = reslice_imgs(deblank(EXPT.FILES.im_files{1}(1,:)),EXPT.(fname),1);
                end
                disp('Checked mask against first image in EXPT.FILES.im_files.  Resliced mask if necessary.');
            catch
                % we don't always have im_files.
            end



            % ---------------------------------------------
            % for fir dxbeta images
            % ---------------------------------------------
        case lower('dxbetas')
            startdir = pwd;
            disp('Are the dx_beta*img files in subject directories below this directory (type 1 or 0):');
            isok = input(startdir);
            if ~isok, startdir = scan_get_files(-1,'*','Select parent dir for Ss dirs with dx_beta*img'); end
            
            fprintf(1,'Getting individual HRF timecourses from dx_beta*img in subject directories.\n');
            fprintf(1,'(And making sure they''re in ascending numeric order).\n');
            EXPT = getfunctnames2(EXPT,'dx*img','FIR.dxbetas','.',startdir);
            fprintf(1,'\n');


            % ---------------------------------------------
            % For group analysis
            % for EXPT.SNPM.P, filenames to put into robfit
            % ---------------------------------------------
        case lower('P')

            disp('Collect images for random effects.')
            if ~isfield(EXPT,'subjects') || isempty(EXPT.subjects), error('You need to create EXPT.subjects first.'); end

            wcard = input('Enter image wildcard (e.g., con*img): ','s');

            %d = dir([EXPT.subjects{1} filesep wcard]);
            % More general: can handle subdirectories prefixed to image
            % wildcard.
            % get names of files with relative paths
            d = filenames([EXPT.subjects{1} filesep wcard]);
            dir_prefix = fileparts(wcard);
            names = cell(length(d), 1);
            
            for i = 1:length(d)
                [dd, ff, ee] = fileparts(d{i});
                names{i} = fullfile(dir_prefix, [ff ee]);
            end


            fprintf('Examining first subject.  Found %3.0f images:\n', length(d));
            fprintf('%s\n', d{:})
            fprintf('Assuming same images exist for other subjects in same dir structure.\n')

            for i = 1:length(d)
                for j = 1:length(EXPT.subjects)

                    if j == 1
                        EXPT.SNPM.P{i} = fullfile(pwd, EXPT.subjects{j}, names{i});
                    else
                        EXPT.SNPM.P{i} = char(EXPT.SNPM.P{i}, fullfile(pwd, EXPT.subjects{j}, names{i}));
                    end
                end
            end

% %             for i = 1:length(d)
% %                 EXPT = getfunctnames2(EXPT,d(i).name,'tmp'); tmp = str2mat(EXPT.tmp{:});
% %                 if ~isempty(varargin)
% %                     EXPT.(varargin{1}).(fname){i} = tmp;
% %                     disp(['Saved images in: EXPT.' varargin{1} '.' fname '{' num2str(i) '}']);
% %                 else
% %                     EXPT.(fname){i} = tmp;
% %                     disp(['Saved images in: EXPT.' fname '{' num2str(i) '}']);
% %                 end
% %             end
% %             EXPT = rmfield(EXPT,'tmp');

            % ---------------------------------------------
            % For group analysis
            % for EXPT.SNPM, get contrast names to put into robfit
            % ---------------------------------------------
        case lower('connames')

            if isempty(varargin), error('You need to specify a substructure for the connames field. e.g., ''connames'',''SNPM'''); end

            disp('Add information about the contrast names.')
            if ~isfield(EXPT,'SNPM'), error('You need to create EXPT.SNPM.  Try EXPT = CreateExpt(''CreateField'',''P'',''SNPM'');');  end
            if ~isfield(EXPT.SNPM,'P') || isempty(EXPT.SNPM.P), error('You need to get image names in EXPT.SNPM.P first. Try EXPT = CreateExpt(''CreateField'',''P'',''SNPM'');'); end

            EXPT.SNPM.connames = [];
            for i = 1:length(EXPT.SNPM.P)
                fprintf(1,'\nImages in this contrast: \n');
                disp(EXPT.SNPM.P{i})
                tmp = input('Enter contrast name for this set (avoid spaces/special chars): ','s');
                EXPT.SNPM.connames = str2mat(EXPT.SNPM.connames,tmp);
            end

            EXPT.SNPM.connames = EXPT.SNPM.connames(2:end,:);
            EXPT.SNPM.connums = 1:length(EXPT.SNPM.P);


        otherwise
            %=======================================================================
            error('Unknown action string')

            %======================================================================

    end


    return




