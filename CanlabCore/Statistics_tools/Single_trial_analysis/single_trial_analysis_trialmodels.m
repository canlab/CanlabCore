function single_trial_analysis_trialmodels(EXPT, submat, trialmodels, varargin)
% single_trial_analysis(EXPT, submat, trialmodels, varargin)
%
% inputs: 
% 'start' followed by start slice
% 
% first run trialmodels = single_trial_setup_models(...), then this
% function

    % Modified by Tor Wager, March 2007

    %-------- SOME PARAMETERS ---------------%
    TR = EXPT.TR;                   %TR
    mask = EXPT.mask;
    nscans = EXPT.SPM.nscan;        %vector containing number of scans in each run
    nruns = length(nscans);         %number of runs
    %----------------------------------------%

    startslice = 1;

    for i = 1:length(varargin)
        if ischar(varargin{i})

            switch varargin{i}
                case {'start', 'startslice','slice'}, startslice = varargin{i+1};
                otherwise
                    disp('Unknown string input.')
            end
        end
    end

    % ---------------------------------------------------------------
    %
    %
    %
    % LOOP THROUGH SUBJECTS
    %
    %
    %
    % ---------------------------------------------------------------

    for s = submat;

        disp(' '); disp(['Processing subject ',num2str(s),'...']);


        % ---------------------------------------------------------------
        % set up images and onsets for this subject
        % ---------------------------------------------------------------

        fprintf('Mask is:\n  %s\n\n', mask);

        disp('Getting image names from EXPT:')
        disp(['Name of EXPT.subjects corresponding to these images: ' EXPT.subjects{s}]);

        % Tor changed this on 1/28/07
        img_names = EXPT.FILES.im_files{s};
        
        %ESTIMATE BETAS
        disp('  ')
        disp('  Starting model fitting.');

        disp('  ')
        
        % Note: tor changed this
        % hipass is only passed in if you want to pre-filter timeseries
        % before estimating trial-level betas or add filtering to a GLM
        % matrix
        % you've already added it to the design matrix, so leave empty

        %%mask = '/Volumes/SCNAlpha/Data_and_Tools/Reward1/scalped_avg152T1_graymatter.img';

        %Pw = whole_brain_filter(img_names,[],TR,[],nruns,[],'glm','trialX',trialX,'bf',bf,'ntrials',ntrials, ...
        %    'dographics',0,'startslice',startslice,'mask',mask);

        Pw = whole_brain_filter(img_names,[],TR,[],nruns,[],'glm', 'trialX', trialmodels, ...
            'dographics',0,'startslice', startslice, 'mask', mask);


    end