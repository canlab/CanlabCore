function nonlin_param_mod_brain(ons, modulator, image_names, SETUP, varargin)
% Nonlinear fits with a parametric modulator on a set of brain images
%
% :Usage:
% ::
%
%     nonlin_param_mod_brain(X, image_names, SETUP, [SETUPional inputs])
%
% :Inputs:
%
%   **ons:**
%        onsets for each condition; one cell per
%        condition, one col. vector per series of onsets
%
%   **modulator:**
%        modulator values for each condition; same
%        format as above
%
%   **image_names:**
%        outcome variable; Images (volume names) for each subject, in
%        string matrix (list of image names); 3-D for now!
%
% :SETUP.(fields):
%
%   **.mask:**
%        name of mask image
%
%   **.preprocX:**
%        flag for whether to HP filter X data
%
%   **.preprocY:**
%        flag for whether to HP filter Y data
%
%   **'nopreproc':**
%        Turn off preproc
%
%   **.TR:**
%        repetition time of volume (image) acquisition
%
%   **.HPlength:**
%        high-pass filter length, in s
%
%   **.scans_per_session:**
%        vector of # volumes in each run, e.g., [128 128 128 128 128]
%
%   **.dummyscans:**
%        indices of images in each run that will be modeled
%        with separate dummy variables
%
%   **.startslice:**
%        starting slice number (to resume analysis)
%
% SETUPional inputs:
% Any of the SETUPions in mediation.m
%
% Also: 'nopreproc' to skip preprocessing (i.e., for trial-level inputs)
%
% ..
%    Tor Wager, May 2008
% ..

    
    % ..
    %    Set up preprocessing
    %    To skip, enter 'nopreproc' as var. arg.
    % ..

    [preprochandle, SETUP] = filter_setup(SETUP, varargin{:});

    % SETUPional: preproc Y?  Should do mostly only if not brain
    % Should not do if using trial-level estimates

    N = size(image_names, 1);  % number of subjects; one cell per subject
    
    % do high-pass filter preprocessing on X, if specified
    if SETUP.preprocX, for i = 1:N, X{i} = preprochandle(X{i}); end, end
    
    if ~SETUP.preprocY 
        preprochandle = []; 
    
    else
        tmp = cell(1, N);
        for i = 1:N, tmp{i} = preprochandle; end
        preprochandle = tmp;
    end

    tr = SETUP.TR;
    
    
    % ---------------------------------------------------------------------
    % Set up mask
    % ---------------------------------------------------------------------
    SETUP.mask_unresampled = SETUP.mask;
    disp('Writing resampled mask.img in current directory.');
    
    scn_map_image(SETUP.mask, image_names(1, :), 'write', 'mask.img');
    SETUP.mask = fullfile(pwd, 'mask.img');
    
    % ---------------------------------------------------------------------
    % Set up analysis
    % ---------------------------------------------------------------------
    
    % THIS stuff is the same for each voxel
    % --------------------------------------
    modulator_centered = modulator - nanmean(modulator);
    xvals = (1:N)';

 
    % This runs the whole thing given a data vector (y)
    fhandle = @(y) nonlin_param_modulator(y, ons, modulator_centered, tr, xvals);

        
    SETUP.names = {'amplitude_mean.img' 'amplitude_by_modulator.img' 'duration_mean.img' 'duration_by_modulator.img' ...
    'intercept.img' 'auc_mean_trial.img' ...
    'auc_by_modulator.img' 'variance.img' };

    SETUP.preprochandle = preprochandle;
    SETUP.fhandle = fhandle;

    SETUP.data.descrip = 'Data after any preprocessing specified (for X and Y)';
    SETUP.data.ons = ons;
    SETUP.data.modulator = modulator;
    SETUP.data.modulator_centered = modulator_centered;
    SETUP.data.image_names = image_names;
    save nonlin_param_mod_SETUP SETUP

    % ---------------------------------------------------------------------
    % Run preprocessing and analysis
    % ---------------------------------------------------------------------
    if ~isfield(SETUP, 'startslice') || isempty(SETUP.startslice), SETUP.startslice = 1; end

    [a, b, c, d, e, f, g, h] = image_eval_function(image_names, fhandle, 'mask', SETUP.mask, 'preprochandle', preprochandle, 'outnames', SETUP.names, 'start', SETUP.startslice);




end    % End Main Function





% Set up preprocessing
function [preprochandle, SETUP] = filter_setup(SETUP, varargin)

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
                    if ~isfield(SETUP, 'TR') || isempty(SETUP.TR), SETUP.TR = NaN; end
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

    SETUP.preproc_any = SETUP.preprocX || SETUP.preprocY

    if ~isfield(SETUP, 'TR') || isempty(SETUP.TR) || isnan(SETUP.TR) 
        error('SETUP.TR must be entered.');
    end
    
    if SETUP.preproc_any && hpflag

        
        error('THIS CODE IS NOT SET UP TO PREPROCESS YET.');
    
        [tmp, I, S] = hpfilter(X(:,1), SETUP.TR, SETUP.HPlength, SETUP.scans_per_session, SETUP.dummyscans); % creates intercept and smoothing matrices

        preprochandle = @(Y) hpfilter(Y,  [], S, SETUP.scans_per_session, I);  % function handle with embedded fixed inputs

    end

end

