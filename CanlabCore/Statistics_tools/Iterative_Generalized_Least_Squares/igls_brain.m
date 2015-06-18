function igls_brain(X, Y, SETUP, varargin)
    % igls_brain(X, Y, SETUP, [SETUPional inputs])
    %
    % Multilevel IGLS/RIGLS on a set of brain images
    %
    % Inputs
    % ------------------------------------------------
    % X                     data matrix of t timepoints x N subjects
    % Y                     outcome variable; Images for each subject, in
    %                       cell array (one subject per cell)
    %
    % SETUP.(fields)
    % .mask                 name of mask image
    % .preprocX             flag for whether to HP filter X data
    % .preprocY             flag for whether to HP filter Y data
    %
    % .TR                   repetition time of volume (image) acquisition
    % .HPlength             high-pass filter length, in s
    % .scans_per_session    vector of # volumes in each run, e.g., [128 128 128 128 128]
    % .dummyscans           indices of images in each run that will be modeled
    %                       with separate dummy variables
    % .startslice           starting slice number (to resume analysis)
    % SETUPional inputs:
    % Any of the SETUPions in mediation.m
    % Also: 'nopreproc' to skip preprocessing (i.e., for trial-level inputs)
    %
    % DESCRIPTION:
    % igls_brain is a "shell" function. It takes IN a set of image names for a voxel-wise analysis.
    % It's job is to specify what the output images should be named, and to
    % specify which function should be run on each voxel.
    % Then it calls a function for generically running analyses across
    % brain voxels and saving maps.  This function is called
    % image_eval_function_multisubj.  
    % So igls_brain just says, run image_eval_function_multisubj.m with
    % THIS igls-specific function and THESE images, and save output images
    % with THESE names
    %
    % For igls analysis, the function is called igls_brain_multilev_wrapper
    % This function, and any other function you want to use a similar
    % "shell" function (like igls_brain) for, must take inputs and outputs in a specific way.
    % Inputs must be defined with fixed parameters, and the only variable
    % parameter (input parameter) is the data vector, Y.  This data vector
    % is replaced with data from each voxel in the brain in the course of
    % making voxel-wise maps. igls_brain sets up a function handle for
    % igls_brain_multilev_wrapper in an appropriate way for igls analysis.
    %
    % You have to pass in:
    % igls_brain(X, Y, SETUP, and then ANY of the igls-specific input
    % arguments that igls.m knows how to understand.)
    %
    % Tor Wager & Martin Lindquist, March 2008
    %
    % Examples:
    %
    %
    %

    

    
    % ---------------------------------------------------------------------
    % Set up preprocessing
    % To skip, enter 'nopreproc' as var. arg.
    % ---------------------------------------------------------------------

    [preprochandle, SETUP] = filter_setup(SETUP, X, varargin{:});

    % SETUPional: preproc Y?  Should do mostly only if not brain
    % Should not do if using trial-level estimates

    N = length(Y);  % number of subjects; one cell per subject
    
    % do high-pass filter preprocessing on X, if specified
    if SETUP.preprocX, for i = 1:N, X{i} = preprochandle(X{i}); end, end
    
    if ~SETUP.preprocY 
        preprochandle = []; 
    
    else
        tmp = cell(1, N);
        for i = 1:N, tmp{i} = preprochandle; end
        preprochandle = tmp;
    end

    % ---------------------------------------------------------------------
    % Set up analysis
    % ---------------------------------------------------------------------

    fhandle = @(Y) igls_brain_multilev_wrapper(X, Y, 'noverbose', varargin{:});


    SETUP.names = {'intercept_b.img' 'slope_b.img' 'intercept_rfxvar_b.img' 'slope_rfxvar_b.img' ...
    'intercept_indiv.img' 'slope_indiv.img' ...
    'sigma.img' 'isconverged.img' 'intercept_t.img' 'slope_t.img'  ...
    'intercept_rfxvar_t.img' 'slope_rfxvar_t.img' 'intercept_p.img' 'slope_p.img'  ...
    'intercept_rfxvar_p.img' 'slope_rfxvar_p.img' , 'intercept_LRTrfxvar_p.img', 'slope_LRTrfxvar_p.img', ...
    'cov_int_b.img' 'cov_int_t.img' 'cov_int_p.img' 'cov_slope_b.img' 'cov_slope_t.img' 'cov_slope_p.img'};

    SETUP.preprochandle = preprochandle;
    SETUP.fhandle = fhandle;

    SETUP.data.descrip = 'Data after any preprocessing specified (for X and Y)';
    SETUP.data.X = X;
    SETUP.data.Y = Y;
    SETUP.data.covariate = [];

    iscovt = strcmp(varargin, 'covariate');
    if any(iscovt)
        iscovt = find(iscovt) + 1;
        disp('2nd-level Covariate found. Maps of covariate stats should be written with meaningful values.')
        SETUP.data.covariate = varargin{iscovt(1)};
        
    else
        disp('2nd-level Covariate not found. Maps of covariate stats will be empty.')
    end
    
    save igls_SETUP SETUP

    % ---------------------------------------------------------------------
    % Run preprocessing and analysis
    % ---------------------------------------------------------------------
    if ~isfield(SETUP, 'startslice') || isempty(SETUP.startslice), SETUP.startslice = 1; end

    image_eval_function_multisubj(Y, fhandle, 'mask', SETUP.mask, 'preprochandle', preprochandle, 'outnames', SETUP.names, 'start', SETUP.startslice);




end    % End Main Function





% Set up preprocessing
function [preprochandle, SETUP] = filter_setup(SETUP, X, varargin)

    preprochandle = [];
    wh_elim = [];
    hpflag = 1;  % only does it if requested, though

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'custompreproc'
                    preprochandle = varargin{i + 1};  % e.g., 'custompreproc', @(data) scale(data) for z=scores;
                
                    hpflag = 0;
                    SETUP.TR = NaN;
                    SETUP.HPlength = [];
                    SETUP.dummyscans = [];
                    wh_elim = i;
                    
                case {'nopreproc'}
                    hpflag = 0;
                    SETUP.preprocY = 0; SETUP.preprocX = 0;
                    SETUP.TR = NaN;
                    SETUP.HPlength = [];
                    SETUP.dummyscans = [];
                    wh_elim = i;
                        
                    % We need to allow mediation SETUPions here, so eliminate this from list and do not error check here.
                    %otherwise, warning(['Unknown input string SETUPion:' varargin{i}]);
            end
        end
    end

    varargin(wh_elim) = [];

    N = fieldnames(SETUP);
    for i = 1:length(N)
        if ~isfield(SETUP, N{i}) || isempty(SETUP.(N{i}))
            switch N{i}
                case {'TR', 'mask', 'scans_per_session', 'preprocY'}
                    error(['Enter SETUP.' N{i}]);

                case 'HPlength'
                    SETUP.(N{i}) = [];

                case 'dummyscans'
                    SETUP.(N{i}) = 1:2;

                otherwise
                    disp('Warning! Unrecognized field in SETUPions structure SETUP.');
            end
        end
    end

    SETUP.preproc_any = SETUP.preprocX || SETUP.preprocY;

    if SETUP.preproc_any && hpflag

        [tmp, I, S] = hpfilter(X(:,1), SETUP.TR, SETUP.HPlength, SETUP.scans_per_session, SETUP.dummyscans); % creates intercept and smoothing matrices

        preprochandle = @(Y) hpfilter(Y,  [], S, SETUP.scans_per_session, I);  % function handle with embedded fixed inputs

    end

end

