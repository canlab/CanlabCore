function trialmodels = single_trial_setup_models(EXPT, onsets, eventdesign, s, varargin)

    % trialmodels = single_trial_setup_models(EXPT, onsets, eventdesign, s, varargin)
    %Constructs a trial-level design matrix, estimates betas and writes them to
    %img files.
    %
    %USAGE: single_trial_analysis(EXPT,onsets,eventdesign,submat);
    %       e.g.,
    %       [onsets2 runlabels eventlabels rundesign eventdesign] = ...
    %       onsets2singletrial(subjects,EXPT.SPM.nscan,'TR',2,eventnames,'noadd');
    %
    %       single_trial_analysis(EXPT,onsets2,eventdesign,1,'gamma3',60);
    %
    %       single_trial_analysis(EXPT,onsets2,eventdesign,2,'gamma3',60, ...
    %       'ratings',ratings,'conditioncontrasts',[-3 -1 1 3], ...
    %       'conditionnames',{'warm' 'low' 'medium' 'hot'}, 'mvmt', ravols_movement_params);
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
    %A useful function in this regard is 'onsets2singletrial',which converts
    %structure-format onsets to cell-format and adjusts the timing. The output of
    %onsets2singletrial can be fed directly into single_trial_analysis.
    %
    % Sam Gershman, 2006
    % Modified by Tor Wager, March 2007

    %-------- SOME PARAMETERS ---------------%
    TR = EXPT.TR;                   %TR
    nscans = EXPT.SPM.nscan;        %vector containing number of scans in each run
    nruns = length(nscans);         %number of runs
    %----------------------------------------%

    %----- SET DEFAULTS ---------------------%
    basistype = 'gamma3';
    hipass = 100;  
    domvmt2 = 0;
    %----------------------------------------%
    

    disp(' '); fprintf('Building design matrix for %3.0f...\n',s);
    disp(['Name of this subject in EXPT.subjects : ' EXPT.subjects{s}]);


    ratings = []; conditioncontrasts = []; conditionnames = []; startslice = 1;

    for i = 1:length(varargin)
        if ischar(varargin{i})

            switch varargin{i}

                case 'ratings', ratings =  varargin{i+1};
                    disp('Found ratings. Will include in model.');
                    
                case {'conditionnames'}, conditionnames = varargin{i+1};
                case {'conditioncontrasts'}, conditioncontrasts = varargin{i+1};

                case {'mvmt','movement'}, mvmt =  varargin{i+1}; disp('Found movement parameters.') 
                        
                case 'mvmt2', domvmt2 = 1; disp('Squared/lagged movement model is ON.'); 
                    
                case {'basis','basisset','bases','basistype'} % name(s) of basis sets to use
                    disp('Found basistype.');
                    basistype = varargin{i+1}; varargin{i+1} = [];
                    
                case 'whichbasis' % enter index in basis sets for each event
                    disp('Found whichbasis.'); 
                    whichbasis = varargin{i+1};
                    
                case {'hipass','hplength','hp'}, hipass = varargin{i+1};
                    disp('Found high-pass filter length.')
                    
                otherwise
                    disp('Unknown string input.')
            end
        end
    end

    


        % ---------------------------------------------------------------
        % set up images and onsets for this subject
        % ---------------------------------------------------------------

        

        ons = onsets{s};    %onsets for subj

        % Tor added this: we now need to enter ntrials to
        % whole_brain_filter
        ntrials = length(ons);

        % ---------------------------------------------------------------
        % CONSTRUCT TRIAL-LEVEL DESIGN MATRIX (models single trials)
        % ---------------------------------------------------------------

        disp('  Constructing design matrix...');

        if ~exist('whichbasis', 'var')
            fprintf('Using the same basis set for all events : %s\n', basistype);

            [trialX, trial_px, bf] = trial_level_beta3('onsets',ons,'output','design', ...
                'rows',sum(nscans),'basistype',basistype,'tr',TR);

            trialX = trialX(:, 1:end-1);    %remove 'experiment' intercept

        else
            fprintf('   Using basis sets indexed by whichbasis for each event\n');
            fprintf('   This will ONLY work if all basis sets have same # parameters.\n');

            trialX = [];

            for i = 1:length(ons)
                
                fprintf(1,'%3.0f ', i);
                
                if i == 1
                    [trialX_thistrial, trial_px, bf] = trial_level_beta3('onsets',ons(i),'output','design', ...
                    'rows',sum(nscans),'basistype',basistype{whichbasis(i)},'tr',TR);
                else
                    [trialX_thistrial] = trial_level_beta3('onsets',ons(i),'output','design', ...
                    'rows',sum(nscans),'basistype',basistype{whichbasis(i)},'tr',TR);
                end

                trialX_thistrial = trialX_thistrial(:,1:end-1); % get rid of intercept

                trialX = [trialX trialX_thistrial];

            end

            fprintf(1, '\n');
        end


        fprintf(1,'   Number of basis functions: %3.0f\n',size(bf,2))
        disp('NOTE: all basis sets for different event types MUST have same number of basis functions!')

        %CONSTRUCT HIGH-PASS FILTERING MATRIX
        % ---------------------------------------------------------------
        disp('  Adding high-pass filter to design matrix...');
        khstr = [];
        for r = 1:nruns
            
            [dummy,KL,KH{r}] = use_spm_filter(TR,nscans(r),'none','specify',hipass);
            if r==1; kh = ['KH{',num2str(r),'}']; else kh = [',KH{',num2str(r),'}']; end
            khstr = strcat(khstr,kh);
            
        end
        KH = eval(['blkdiag(',khstr,')']);


        % Add movement param stuff here
        % ---------------------------------------------------------------
        MV = [];
        if exist('mvmt', 'var')

            disp('  Adding movement parameters...');
            
            if size(mvmt, 1) ~= size(trialX, 1)
                error('mvmt params and trial design are different lengths.');
            end

            mvmt = scale(mvmt);
            
            if domvmt2
                disp('   Adding lagged and squared movement parameters...');
                m2 = [zeros(1, 6); mvmt(1:end-1,:)];
                %gradient(mvmt);
                MV = [mvmt mvmt .^ 2 m2 m2 .^ 2];

            else
                MV = mvmt;
            end

        end


        % Add intercepts and dummy covs for the first 2 trials
        % ---------------------------------------------------------------
        disp('  Adding intercepts and dummy covs for first 1 trials...');
        whremove = 1;  % could be 1:2, 1:3, etc.
        I = intercept_model(nscans, 1, whremove);
        
        
        % Put all regressors together
        % ---------------------------------------------------------------  
        trialX_of_interest = trialX;
        
        trialX = [trialX KH MV I];
        trialpx = pinv(trialX);

        % save
        % ---------------------------------------------------------------
        disp(['  Saving to: ',EXPT.subjects{s},'_design.mat']);
        save([EXPT.subjects{s},'_design'], 'trialX', 'trialpx', 'trialX_of_interest');




        % ---------------------------------------------------------------
        % Get info for SUBJECT-LEVEL DESIGN MATRIX (models condition
        % effects across trials)
        % ---------------------------------------------------------------
        
        % Each row is a trial, each col. is a task condition
        % indicator matrix of 1/0, 1's in col. k where trials belong to condition k
        
        % Build subject-level design matrix, for analysis of condition
        % effects across trials
        % No intercept: interpret change from abs. baseline of zero
        
        subjlevelX = eventdesign{s};

        if size(subjlevelX, 1) ~= ntrials
            error('eventdesign is wrong size for ntrials.  Error in either onsets2singletrial or ons.');
        end
        
        subjlevelpx = pinv(subjlevelX);

        % contrasts should be in column vectors
        if isrow(conditioncontrasts), conditioncontrasts = conditioncontrasts'; end

        % stuff needed for whole_brain_filter trial-level models
        % ---------------------------------------------------------------
        
        trialmodels.trialX = trialX;
        trialmodels.bf = bf;
        trialmodels.ntrials = ntrials;
        trialmodels.trialpx = trialpx;
        trialmodels.subjlevelX = subjlevelX;
        trialmodels.subjlevelpx = subjlevelpx;
        trialmodels.conditioncontrasts = conditioncontrasts;
        trialmodels.conditionnames = conditionnames;
        
        if exist('ratings','var')
            disp('Found ratings and will include in trialmodels var, which is used in model fitting.');
            disp('So: will analyze by ratings in model fitting.');
            
            trialmodels.ratings = ratings{s}; 
        end
        
        % stuff to save even tho we don't need it in fitting.
        trialmodels.ons = ons;
        trialmodels.basistype = basistype;
        trialmodels.ntrials = ntrials;
  

        save([EXPT.subjects{s},'_design'],'-append','trialmodels');
        
        % tor added this to make image; 
        figure; imagesc(trialX), colormap gray; drawnow


disp('Done setting up models.');


    end