%Constructs a trial-level design matrix, estimates betas and writes them to
%img files.
%
%USAGE: single_trial_analysis(EXPT, onsets, eventdesign, submat);
%       e.g.,
%       [onsets2 runlabels eventlabels rundesign eventdesign] = ...
%       onsets2singletrial(subjects, EXPT.SPM.nscan, 'TR', 2, eventnames, 'noadd');
%
%       single_trial_analysis(EXPT, onsets2, eventdesign, 1, 'gamma3', 60);
%
%       single_trial_analysis(EXPT, onsets2, eventdesign, 2, 'gamma3', 60, ...
%       'ratings', ratings, 'conditioncontrasts', [-3 -1 1 3], ...
%       'conditionnames', {'warm' 'low' 'medium' 'hot'}, 'mvmt', ravols_movement_params);
%
%INPUTS:
%   EXPT - EXPT object
%       Required fields:
%       subjects
%       TR
%       mask
%       SPM.nscan
%       FILES.im_files
%
%   onsets - 'cell' format onsets (see below)
%       to convert from structure format, see onsets2singletrial
%
%   eventdesign - stick (delta) function indicating condition onsets;
%                 returned as output of onsets2singletrial
%
%   submat - subject (or vector of subject) indices
%   basistype (optional) - basis type (default: 'gamma3')
%   hipass (optional) - high-pass filter (default: 60)
%   startslice (optional) - starting slice
%
% Optional inputs:
% conditioncontrasts: followed by contrasts across expt'l conditions
% which are coded in eventdesign.  e.g., [-3 -1 1 3]';
%
% conditionnames: names for each condition in eventdesign: {'warm' 'low' 'medium' 'hot'};
%
% ratings : followed by behavioral data values
% in same format as onsets
% See onsets2singletrial for ways to do this.
%
% stimlength: the length of the pain stimulus (in seconds) for use with
%   basis type "painX"
%
%OUTPUTS:
%   -saves a mat file with the suffix 'design' that contains the design mtx
%
%see trial_level_beta3.m for more information
%
%A NOTE ABOUT ONSETS:
%The SCAN lab uses two different formats for its onsets. The
%'structure' format contains a different field for each subject, and a
%different sub-field for each run, with sub-fields for each condition that
%contain the onsets. The 'cell' format is a cell vector, where each cell
%contains 1 subject; within each subject-cell is another cell vector, where
% each cell contains onsets for a particular run/condition. If, for example, you have 2
%conditions and 2 runs, the cell array will be organized as follows:
%The 1st cell will contain the onsets of all events from condition 1 in
%run 1, the 2nd cell will contain the onsets of events from condition 2 in
%run 1, the 3rd cell will contain onsets of events from condition 1 in run
%2, etc. Specifically for single-trial analysis, the onsets must be in cell
%format (and in TRs), with the added requirement that the onsets be adjusted
%so that the first onset of each run is equal to the onset of the last scan
%(*not event*) of the previous run. For example, if you have 100 scans in
%each run, the first onset of the 1st run should be 0, the first onset of
%the 2nd run should be 100, etc.
%A useful function in this regard is 'onsets2singletrial', which converts
%structure-format onsets to cell-format and adjusts the timing. The output of
%onsets2singletrial can be fed directly into single_trial_analysis.
%
% Sam Gershman, 2006
% Modified by Tor Wager, March 2007

function single_trial_analysis(EXPT, onsets, eventdesign, submat, basistype, hipass, varargin)
    covs = [];
    
    TR = EXPT.TR;                   
    mask = EXPT.mask;
    nscans = EXPT.SPM.nscan;        %vector containing number of scans in each run
    nruns = length(nscans);         %number of runs
    if(isfield(EXPT, 'cov') && ~isempty(EXPT.cov))
        covs = EXPT.cov;
    end

    %----- SET DEFAULTS ---------------------%
    if nargin < 5, basistype = 'gamma3'; hipass = 60; end
    if nargin < 6, hipass = 60;  end


    ratings = []; conditioncontrasts = []; conditionnames = []; startslice = 1;
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case {'start', 'startslice'}
                    startslice = varargin{i+1};
                case 'ratings'
                    ratings =  varargin{i+1};
                case 'conditionnames'
                    conditionnames = varargin{i+1};
                case 'conditioncontrasts'
                    conditioncontrasts = varargin{i+1};
                case {'mvmt', 'movement'}
                    mvmt = varargin{i+1};
                case 'spikes'
                     spikecovs = varargin{i+1};
                case {'covs', 'nuisancecovs'}
                    if(~isempty(covs) && covs ~= varargin{i+1})
                        error('The covariates in EXPT.cov and the covariates passed in do not match.');
                    else
                        covs = varargin{i+1};
                    end
                case 'stimlength'
                    stimlength = varargin{i+1};

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

    for s = submat
        disp(' ');
        disp(['Processing subject ', num2str(s), '...']);


        % ---------------------------------------------------------------
        % set up images and onsets for this subject
        % ---------------------------------------------------------------

        fprintf('Mask is:\n  %s\n\n', mask);

        disp('Getting image names from EXPT:')
        disp(['Name of EXPT.subjects corresponding to these images: ' EXPT.subjects{s}]);

        % Tor changed this on 1/28/07
        img_names = EXPT.FILES.im_files{s};
        % Before 1/28/07 THIS WAS: img_names = EXPT.FILES.im_files{1};     %vector of image names for subj

        ons = onsets{s};    %onsets for subj

        % Tor added this: we now need to enter ntrials to
        % whole_brain_filter
        ntrials = length(ons);

        % ---------------------------------------------------------------
        % CONSTRUCT TRIAL-LEVEL DESIGN MATRIX (models single trials)
        % ---------------------------------------------------------------
        disp('  Constructing design matrix...');
        if strcmp(basistype,'painX')
            if ~exist('stimlength','var')
                error('basis type ''painX'' requires specification of stimlength')
            end
            [trialX, trial_px, bf] = trial_level_beta3('onsets', ons, 'output', 'design', 'rows', sum(nscans), 'basistype', basistype, 'tr', TR, 'stimlength', stimlength);
        else
            [trialX, trial_px, bf] = trial_level_beta3('onsets', ons, 'output', 'design', 'rows', sum(nscans), 'basistype', basistype, 'tr', TR);
        end
        trialX = trialX(1:end, 1:end-1);    %remove 'experiment' intercept

        fprintf('Number of basis functions: %3.0f\n', size(bf, 2))

        %CONSTRUCT HIGH-PASS FILTERING MATRIX
        % ---------------------------------------------------------------
        disp('  Adding high-pass filter to design matrix...');
        clear KH
        for r = 1:nruns
            [dummy, KL, KH{r}] = use_spm_filter(TR, nscans(r), 'none', 'specify', hipass);
        end
        KH = blkdiag(KH{:});


        % Add movement param, spike, and other covariate stuff here
        % ---------------------------------------------------------------
        % check length
        disp('  Checking images and covariates...');
        allimgs = expand_4d_filenames(img_names);
        
        nimgs = size(allimgs, 1);
        fprintf('\nFound %3.0f images\n', nimgs);
        
        if nimgs ~= sum(nscans)
            error('Num images does not match number of scans entered for each run.');
        end
        
        MV = [];
        if exist('mvmt', 'var')
            fprintf('Movement parameters found.  Adding squared, shifted, squared gradient params.\n');
            if size(mvmt, 1) ~= nimgs
                error('mvmt params and images are different lengths.');
            end

            mvmt = scale(mvmt);
            m2 = [zeros(1, 6); mvmt(1:end-1,:)];
            %gradient(mvmt);
            MV = [mvmt mvmt.^2 m2 m2.^2];
        else
            fprintf('Movement parameters not entered.\n');
        end

        if exist('spikecovs', 'var')
            fprintf('Spikes/additional covariates found.\n');
            if size(spikecovs, 1) ~= nimgs
                error('spike covs and images are different lengths.');
            end
            
            MV = [MV spikecovs];
        else
            fprintf('Spikes/additional covariates not entered.\n');
        end

        % Add intercepts and dummy covs for the first 2 trials
        % ---------------------------------------------------------------
        disp('  Adding intercepts and dummy covs for first 1 trials...');
        whremove = 1;  % could be 1:2, 1:3, etc.
        I = intercept_model(nscans, 1, whremove);


        % Put all regressors together
        % ---------------------------------------------------------------
        trialX_of_interest = trialX; %#ok;

        trialX = [trialX KH MV scale(covs) I];
        trialpx = pinv(trialX);

        % save
        % ---------------------------------------------------------------
        disp(['  Saving to: ', EXPT.subjects{s}, '_design.mat']);
        save([EXPT.subjects{s}, '_design'], 'trialX', 'trialpx', 'trialX_of_interest');




        % ---------------------------------------------------------------
        % Get info for SUBJECT-LEVEL DESIGN MATRIX (models condition
        % effects across trials)
        % ---------------------------------------------------------------

        % Each row is a trial, each col. is a task condition
        % indicator matrix of 1/0, 1's in col. k where trials belong to condition k

        % Build subject-level design matrix, for analysis of condition
        % effects across trials
        % No intercept: interpret change from abs. baseline of zero

        evt = eventdesign{s};
        subjlevelX = evt;

        subjlevelpx = pinv(subjlevelX);

        % contrasts should be in column vectors
        if isrow(conditioncontrasts)
            conditioncontrasts = conditioncontrasts';
        end

        % stuff needed for whole_brain_filter trial-level models
        % ---------------------------------------------------------------
        trialmodels.trialX = trialX;
        trialmodels.bf = bf;
        trialmodels.basistype = basistype;
        trialmodels.ntrials = ntrials;
        trialmodels.trialpx = trialpx;
        trialmodels.subjlevelX = subjlevelX;
        trialmodels.subjlevelpx = subjlevelpx;
        trialmodels.conditioncontrasts = conditioncontrasts;
        trialmodels.conditionnames = conditionnames;
        if ~isempty(ratings)
            trialmodels.ratings = ratings{s};
        end

        save([EXPT.subjects{s} '_design'], '-append', 'trialmodels');

        %Canonical model; do we need this?
        % %         evt = eventdesign{s};
        % %         nconds = size(evt, 2);
        % %         delta = zeros(size(trialX, 1), nconds);
        % %         delta(ons+1,:) = evt;
        % %         canon = getPredictors(delta, spm_hrf(TR) ./ max(spm_hrf(TR)));    %canon is convolved predictors, collapsing across runs

        % tor added this to make image; are all Ss onsets different?
        create_figure('Design matrix'); set(gca,'YDir','Reverse');
        imagesc(trialX); axis auto, axis tight
        colormap gray; 
        drawnow

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

        %Pw = whole_brain_filter(img_names, [], TR, [], nruns, [], 'glm', 'trialX', trialX, 'bf', bf, 'ntrials', ntrials, ...
        %    'dographics', 0, 'startslice', startslice, 'mask', mask);

        Pw = whole_brain_filter(img_names, [], TR, [], nruns, [], 'glm', 'trialX', trialmodels, ...
            'dographics', 0, 'startslice', startslice, 'mask', mask);
    end
end
