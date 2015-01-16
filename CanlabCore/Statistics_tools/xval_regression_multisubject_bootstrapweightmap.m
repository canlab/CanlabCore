function STATS = xval_regression_multisubject_bootstrapweightmap(fit_method, Y, X, volInfo, varargin)
    % STATS = xval_regression_multisubject_bootstrapweightmap(fit_method, Y, X, volInfo, varargin)
    %
    % CROSS-VALIDATED (JACKKNIFE) REGRESSION
    % Leave-one observation out, predict outcomes for each missing holdout_set.
    % THIS IS LIKE XVAL_REGRESSION_MULTISUBJECT, but IS SPECIFIC FOR
    % CREATING P-VALUES FOR BOOSTRAPPED VOXEL WEIGHTS. So it requires an
    % additional input, volInfo (see below).
    %
    % Y = outcome, holdout_set x 1
    % X = predictors, holdout_set x variables
    % fit_method = 'lasso' 'ols' 'ridge' or 'robust'
    % optional inputs:
    % 'pca', PCA data reduction first
    % 'ndims', dims to save from PCA
    % 'variable', retain max dims for each subject
    % case {'cov', 'covs'}, cov = varargin{i+1};  NOTE:covariates are added to training data if entered.
    % {'lassopath', 'ridgek', 'regparams'}, regparams = varargin{i+1};
    % {'noverb', 'noverbose'}, verbose = 0;
    % {'lowverb', 'loverb', 'lowverbose'}, verboseL = 1; verbose = 0;
    % 'nested_choose_ndims', dochoose_ndims = 1;
    %
    % volInfo: volume info structure with locations of voxels in X
    %            as created by iimg_read_img
    %
    % Single-level analysis: enter a single cell with Y and X data.
    % rows are observations. (e.g., subjects)
    %
    % Multi-level analysis: enter a cell per subject with Y and X data.
    % rows are observations within-subjects and should be independent for valid
    % cross-validation.
    %
    % STATS = xval_regression_multisubject('lasso', pain_scores, data, 'pca',
    % 'ndims', 'variable');
    % SEE THE WIKI FOR EXAMPLES, ETC.
    % See also xval_lasso_brain.m
    % See also xval_regression_bootstrapweightmap_brainplots(STATS)
    %
    % Tor Wager, 3/17/09

    % Set defaults
    % -----------------------------------------------------------
    nested_setdefaults();
    verbose;

    % Variable dimensions: retain max possible for each subject
    % If choose ndims is selected, then pick the number here
    % Based on size of input data only (not optimized; full dims)
    % -----------------------------------------------------------
    %nested_choose_ndims();

    [STATS.INPUTS.pcsquash, STATS.INPUTS.num_dims, STATS.INPUTS.Y, STATS.INPUTS.holdout_method] = deal(pcsquash, num_dims, Y, holdout_method);
    %STATS.INPUTS.X = X;

    %include = 1:N;  % which subjects to include

    for s = 1:N
        % ===========================================
        % DATASET (SUBJECT) LOOP
        % This function will run separately for each dataset ("subject" in a
        % multi-level analysis) in the input cell arrays
        % With one dataset, it runs cross-validation for that dataset only
        % ===========================================

        nested_prepdata(); % remove NaNs and initialize fit variable (output)
        nanvox; % we will need to pass this into another inline later

        % Return holdout_set{} defining folds and holdout set for each fold
%         holdout_set = nested_select_holdout_set();
%         STATS.INPUTS.holdout_set{s} = holdout_set;

        train_y = Y{s}; % all data, here -- not leaving out holdout set

        if ~isempty(cov_val)
            train_dat = [X{s} cov_val{s}];
        else
            train_dat = [X{s}];
        end

        test_dat = train_dat; % no holdout; this is to get voxel weights

        % Dim reduction
        % --------------------------------------------------------------------------------
        if pcsquash
            [v, train_dat, test_dat] = do_pcsquash(train_y, train_dat, test_dat, num_dims(s), doplssquash);
            %STATS.VOXWEIGHTS.eigenvectors = v;
        end

        % Fit
        % --------------------------------------------------------------------------------
        % subjbetas{s} is a list of voxel weights for one subject, typically
        nvox = size(X{s}, 2);  % original voxels, not including covs or NaN voxels (will add in NaNs later)
        [STATS.VOXWEIGHTS.betas{s}, STATS.VOXWEIGHTS.vox_weights{s}] = do_fit(fit_method, train_y, train_dat, pcsquash, v, nvox, regparams);

        % --------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------

        % Bootstrap the above

        % --------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------

        nbootsamples = 5000; % 1000 does not give p < .001.  2001+ does.
        
        bootsam = setup_boot_samples(train_y, nbootsamples);

        STATS.INPUTS.nbootsamples = nbootsamples;
        STATS.VOXWEIGHTS.bootbetas = cell(1, nbootsamples);
        STATS.VOXWEIGHTS.boot_vox_weights = zeros(length(STATS.VOXWEIGHTS.vox_weights{s}), nbootsamples, 'single');
        STATS.VOXWEIGHTS.volInfo = volInfo;
        
        fprintf('Getting p-values for voxel weights: Running %3.0f bootstrap samples: 000', nbootsamples);

        poolstate = matlabpool('size');
        if poolstate  
            fprintf('Parallel mode\n');
            nested_bootstrap_loop_parallel
        else
            nested_bootstrap_loop
        end
        
        testvalue = 0;
        % % i = 1
        bstat = STATS.VOXWEIGHTS.boot_vox_weights'; %(i, :)

        pct_lowertail = sum(bstat < testvalue) ./ nbootsamples;
        pct_uppertail = sum(bstat > testvalue) ./ nbootsamples;

        % lower tail is smaller for positive effect. is_lowertail == 1 effects
        % should be positive.
        is_lowertail = pct_lowertail < pct_uppertail;

        % make sure pct is not 1 or 0 due to limited bootstrap samples
        pct_lowertail = max(pct_lowertail, 1 ./ nbootsamples);
        pct_uppertail = min(pct_uppertail, 1 - 1./nbootsamples);

        pct_uppertail = max(pct_uppertail, 1 ./ nbootsamples);
        pct_lowertail = min(pct_lowertail, 1 - 1./nbootsamples);

        z_lowertail = norminv(pct_lowertail);  % neg z for pos effect, and vice versa
        z_uppertail = norminv(pct_uppertail);

        z = z_uppertail;
        z(is_lowertail) = -z_lowertail(is_lowertail); % sign of this is irrelevant because we take min  (p, 1-p) below

        p = normcdf(z);
        p = [p; 1 - p];
        p = 2 .* min(p);

        % pvals will be min for out-of-mask voxels
        p(STATS.VOXWEIGHTS.vox_weights{s}' == 0) = 1;
        p_for_fdr = p(STATS.VOXWEIGHTS.vox_weights{s}' ~= 0);
        
        STATS.VOXWEIGHTS.pval_2tail = p;
        STATS.VOXWEIGHTS.pval_2tail_descrip = 'Bootstrapped p-values for voxel weights';

        STATS.VOXWEIGHTS.pthr_fdr = FDR(p_for_fdr, .05);
        if isempty(STATS.VOXWEIGHTS.pthr_fdr), STATS.VOXWEIGHTS.pthr_fdr = 0; end
        STATS.VOXWEIGHTS.sigfdr = p < STATS.VOXWEIGHTS.pthr_fdr;

        % ****edit this to make output images
        nested_output_metrics_and_plots();

    end % dataset/subject loop


    STATS.Y_orig = Y_orig;
    if ~isempty(cov_val)
        STATS.note = 'covariates were added to predictive model'; %, if covs %removed in Y_orig stored here, if covs were entered';
    else
        STATS.note = 'no covariates entered';
    end
% % 
% %     STATS.devs_from_mean_only_model = devs_from_mean_only_model;
% %     STATS.devs_from_full_model = devs_from_full_model;
% %     STATS.var_null = var_null;
% %     STATS.var_full = var_full;
% % 
% %     % Null model leave-one-out prediction error
% %     % By predicting based on the mean of OTHER subjects, we've increased
% %     % the variance
% %     STATS.pred_err_null = var_null .^ .5;
% % 
% %     STATS.pred_err = pred_err;
% %     STATS.pred_err_descrip = 'Apparent error/loss: Root mean squared deviation from observed outcomes';
% %     STATS.var_reduction = rsq;
% %     STATS.var_reduction_descrip = 'Percent reduction in outcome variance due to model; negative values indicate added variance.';
% % 
% %     STATS.subjfit = subjfit;
% %     STATS.subjbetas = subjbetas;
% %     STATS.r_each_subject = rr;
% %     STATS.r_squared = STATS.r_each_subject .^ 2;
% %     STATS.r_each_subject_note = 'r value: correlation between predicted and outcome values';
% %     STATS.r_each_subject_note2 = 'may not be very interpretable, because null model r = -1.0';
% % 
% %     if verbose
% %         disp(['Mean correlation is: ' num2str(mean(rr))])
% %         disp(['Mean proportion of original variance explained is: ' num2str(mean(rsq))])
% %         disp('The number above is based on reduction of original variance;');
% %         disp('it can be negative if predictors are not helpful because they add noise, increasing the overall variance!')
% %         disp(' ')
% %         fprintf('Null model pred. error is %3.2f, and full model is %3.2f\n', STATS.pred_err_null(1), STATS.pred_err(1));
% %         disp(' ')
% %     end



    % =======================================================================
    % =======================================================================
    %
    %Inline (nested) functions
    %
    % =======================================================================
    % =======================================================================

    function nested_setdefaults
        pcsquash = 0;
        doplssquash = 0;
        num_dims = 2;
        cov_val = [];
        verbose = 1;
        verboseL = 0;
        lassopath = '/Users/tor/Documents/matlab_code_external/machine_learning/lasso_rocha/lasso';
        doplot = 1;
        dochoose_ndims = 0;
        regparams = [];
        dochoose_regparams = 0;
        holdout_method = 'loo';

        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}
                    % Dimension reduction
                    case {'pca', 'pcsquash'}, pcsquash = 1;
                    case {'pls', 'plssquash'}, pcsquash = 1; doplssquash = 1;
                    case {'num_dims', 'ndims'}, num_dims = varargin{i+1}; varargin{i + 1} = [];
                    case 'nested_choose_ndims', dochoose_ndims = 1;

                        % Covariates
                    case {'cov_val', 'covs', 'cov'}, cov_val = varargin{i+1};

                        % Shrinkage/regularization parameters
                    case {'lassopath', 'ridgek', 'regparams'}, regparams = varargin{i+1}; varargin{i + 1} = [];
                    case {'choose_regparams', 'dochoose_regparams', 'optimize_regularization'}, dochoose_regparams = 1;

                        % Holdout set selection
                    case {'holdout_method'}, holdout_method = varargin{i+1}; varargin{i + 1} = [];

                        % Output control
                    case {'noverb', 'noverbose'}, verbose = 0;
                    case {'lowverb', 'loverb', 'lowverbose'}, verboseL = 1; verbose = 0;
                    case 'verbose' % default

                    otherwise, warning(['Unknown input string option:' varargin{i}]);
                end
            end
        end

        if verbose
            fprintf('xval_regression_multisubject\nVerbose mode (enter ''lowverbose'' to minimize or ''noverbose'' to turn off verbose output)\n')
        end

        N = length(Y); % number of subjects/datasets
        [subjbetas, subjfit] = deal(cell(1, N));

        fit = zeros(length(Y{1}));

        if strcmp(fit_method, 'lasso'), addpath(lassopath), end

    end % defaults setting function


    % =======================================================================

    function nested_prepdata
        if verbose
            fprintf('Dataset %3.0f\n> -------------------------------\n ', s);
        elseif verboseL
            fprintf('%3.0f ', s);
        end

        % remove NaNs
        % ---------------------------------------------------------------------
        % X variables are done first, on the presumption that variables
        % (voxels) are many and observations are few
        nanvox{s} = any(isnan(X{s}));
        if all(nanvox{s}), error('All X variables appeared to have NaN values for one or more observations.'); end
        if verbose && sum(nanvox{s}), fprintf('Removed %3.0f X variables with NaNs\n', sum(nanvox{s})); end
        X{s}(:, nanvox{s}) = [];

        if isempty(cov_val)
            [wasnan{s}, Y{s}, X{s}] = nanremove(Y{s}, X{s});
        else
            [wasnan{s}, Y{s}, X{s}, cov_val{s}] = nanremove(Y{s}, X{s}, cov_val{s});
        end
        if all(wasnan{s}), error('All observations appeared to have NaN values for one or more variables.'); end
        if verbose && sum(wasnan{s}), fprintf('Removed %3.0f observations with NaNs\n ', sum(wasnan{s})); end

        Y_orig{s} = Y{s};

        fit = NaN * zeros(size(Y{s}));
        tic

        % ---------------------------------------------------------------------
        % initialize optional things
        v = [];

        if verbose || verboseL, fprintf('Fold %3.0f', 0); end

    end


    % =======================================================================

    function nested_choose_ndims()

        if isstr(num_dims) && strcmp(num_dims, 'variable')
            if pcsquash == 0, warning('xval:ConflictingInputs', 'PC squash is off, so num_dims input will not be used.'); end

            if verbose, fprintf('Variable number of dimensions: choosing: '); end

            clear num_dims
            for i = 1:N
                if ~isempty(cov_val)
                    num_dims(i) = min(size(X{i}, 1) - 2, size(X{i}, 2) + size(cov_val{i}, 2) - 2);
                else
                    num_dims(i) = min(size(X{i}, 1) - 2, size(X{i}, 2) - 2);
                end
                if verbose, fprintf('%03d ', num_dims(i)); end
            end
            if verbose, fprintf('\n'); end
        end

        if length(num_dims) == 1 && N > 1, num_dims = repmat(num_dims, N, 1); end

        if dochoose_ndims
            disp('Inner cross-validation; this could take a long time')
        end

    end

    % =======================================================================
% % % 
% % %     function holdout_set = nested_select_holdout_set
% % %         % purpose: return holdout_set variable
% % % 
% % %         nobs_s = length(Y{s});
% % % 
% % %         switch lower(holdout_method)
% % %             case 'loo'
% % %                 holdout_set = cell(1, nobs_s);
% % %                 for i = 1:nobs_s, holdout_set{i} = i; end
% % % 
% % %             case 'l2o'
% % % 
% % %             case 'balanced4'
% % %                 nfolds = 4;
% % %                 holdout_set = cell(1, nfolds);
% % %                 [ys, indx] = sort(Y{s});
% % %                 for k = 1:nfolds
% % % 
% % %                     holdout_set{k} = indx(k:nfolds:end);
% % %                     if isempty(holdout_set{k}), error('Holdout set construction error: Dataset too small?'); end
% % % 
% % %                 end
% % % 
% % %             otherwise error('Unknown holdout method.  See help.');
% % %         end
% % %     end
% % % 
% % % 
% % %     % =======================================================================
% % % 
% % %     function nested_select_training_test_data
% % %         % Select training/test data
% % %         % -----------------------------
% % %         train_y = Y{s};
% % %         train_y(holdout_set{wh_fold}) = [];
% % % 
% % %         if ~isempty(cov_val)
% % %             train_dat = [X{s} cov_val{s}];
% % %         else
% % %             train_dat = [X{s}];
% % %         end
% % % 
% % %         test_dat = train_dat(holdout_set{wh_fold}, :);
% % %         train_dat(holdout_set{wh_fold}, :) = [];              % leave out the missing observation(s)
% % % 
% % %     end % data selection function

    % =======================================================================

    function nested_bootstrap_loop
        
        tic
        
        for booti = 1:nbootsamples

            fprintf('\b\b\b%3.0f', booti);

            train_y = Y{s}(bootsam(:, booti)); % all data, here -- not leaving out holdout set

            if ~isempty(cov_val)
                train_dat = [X{s}(bootsam(:, booti), :) cov_val{s}(bootsam(:, booti), :)];
            else
                train_dat = [X{s}(bootsam(:, booti), :)];
            end

            test_dat = train_dat; % no holdout; this is to get voxel weights

            if pcsquash
                my_ndims = min(length(unique(bootsam(:, booti))) - 2, num_dims(s));
                [v, train_dat, test_dat] = do_pcsquash(train_y, train_dat, test_dat, my_ndims, doplssquash);
                %STATS.VOXWEIGHTS.eigenvectors = v;
            end

            nvox = size(X{s}, 2);  % original voxels, not including covs or NaN voxels (will add in NaNs later)
            [STATS.VOXWEIGHTS.bootbetas{booti}, STATS.VOXWEIGHTS.boot_vox_weights(:, booti)] = do_fit(fit_method, train_y, train_dat, pcsquash, v, nvox, regparams);

        end % end bootstrap

        fprintf('\n')
        toc
        
    end
    % =======================================================================
    function nested_bootstrap_loop_parallel
        
        % slice vars
%        clear train_y train_dat test_dat
%         for booti = 1:nbootsamples
%             
%             
%             
%             % do this here just so we can clear bootsam
%             if pcsquash
%                 my_ndims{booti} = min(length(unique(bootsam(:, booti))) - 2, num_dims(s));
%             end
%             
%             
%         end
        
%         clear bootsam
%        test_dat = train_dat; % no holdout; this is to get voxel weights
        nvox = size(X{s}, 2);  % original voxels, not including covs or NaN voxels (will add in NaNs later)

        [train_dat test_dat train_y v my_ndims] = deal(cell(1, nbootsamples));
        
        % PARALLEL BOOTSTRAP LOOP
        tic
        fprintf('Done: ')
        parfor booti = 1:nbootsamples

            if mod(booti, 10) == 0, fprintf(' %3.0f', booti); end
            drawnow
            
            train_y{booti} = Y{s}(bootsam(:, booti)); % all data, here -- not leaving out holdout set
            
            if ~isempty(cov_val)
                train_dat{booti} = [X{s}(bootsam(:, booti), :) cov_val{s}(bootsam(:, booti), :)];
            else
                train_dat{booti} = [X{s}(bootsam(:, booti), :)];
            end
            
            test_dat{booti} = train_dat{booti}; % no holdout; this is to get voxel weights
            
            if pcsquash
                my_ndims{booti} = min(length(unique(bootsam(:, booti))) - 2, num_dims(s));
                [v{booti}, train_dat{booti}, test_dat{booti}] = do_pcsquash(train_y{booti}, train_dat{booti}, test_dat{booti}, my_ndims{booti}, doplssquash);
            end

            [bootbetas{booti}, boot_vox_weights(:, booti)] = do_fit(fit_method, train_y{booti}, train_dat{booti}, pcsquash, v{booti}, nvox, regparams);

            train_y{booti} = [];
            train_dat{booti} = [];
            test_dat{booti} = [];
            v{booti} = [];
            
        end % end bootstrap

        fprintf('\n')
        toc

        STATS.VOXWEIGHTS.bootbetas = bootbetas;
        STATS.VOXWEIGHTS.boot_vox_weights = boot_vox_weights;
        
    end
    
    
    
        % =======================================================================

        
    function nested_output_metrics_and_plots
        
        % Print summary table
        % ----------------------------------------------
        diary bootstrap_weightmap_output.txt
        
        thresh = [Inf .001 .002 .005 .01 .05];
        sigvox = sum(STATS.VOXWEIGHTS.sigfdr);
        threshnames = {'FDR' '001' '002' '005' '01' '05'};
        
        thresh(1) = STATS.VOXWEIGHTS.pthr_fdr;
        
        for ii = 2:length(thresh)
            sigvox(ii) = sum(STATS.VOXWEIGHTS.pval_2tail < thresh(ii));
        end
        
        fprintf('In-Mask Space is %3.0f\n', volInfo.n_inmask)
        fprintf('FDR threshold = %3.6f, sig vox = %3.0f\n', STATS.VOXWEIGHTS.pthr_fdr, sum(STATS.VOXWEIGHTS.sigfdr))
        for ii = 2:length(thresh)
            fprintf('p < %3.4f, sig vox = %3.0f\n', thresh(ii), sigvox(ii))
        end
        
        diary off
        
        % Set Colors
        % ----------------------------------------------
        reversecolors = 0;
        if reversecolors, revstr = '_revcolors'; else revstr = ''; end
        
        if reversecolors
            negcm = colormap_tor([.8 0 0], [1 1 0], [1 1 0]);
            poscm = colormap_tor([.4 .8 .6], [0 0 1], [0 0 1]);
        else
            poscm = colormap_tor([.8 0 0], [1 1 0], [1 1 0]);
            negcm = colormap_tor([.4 .8 .6], [0 0 1], [0 0 1]);
        end
        
        % Save weight image
        % ----------------------------------------------
        w = STATS.VOXWEIGHTS.vox_weights{1};
        iimg_reconstruct_vols(w, volInfo, 'outname', 'voxelweights.img');
        
        iimg_reconstruct_vols(p', volInfo, 'outname', 'voxelweight_pvalues.img');
        
        % ***save legend, summarize pos and neg weights***
        
        % -------------------------------------------------------------
        % Save weight images and display montages for each threshold
        % -------------------------------------------------------------
        
        % NEW MONTAGE SETUP
        % -------------------------------
        obj = fmridisplay;            % Create object with canonical underlay
        obj = montage(obj, 'onerow');  % Show axial montage of underlay
        %obj = montage(obj, 'axial', 'slice_range', [-40 55], 'onerow', 'spacing', 8);
        enlarge_axes(gcf, .95);
        obj = montage(obj, 'saggital', 'slice_range', [-6 6], 'onerow');

        for i = 1:length(thresh)
            
            if ~sigvox(i), continue, end  % skip if no vox
            
            outbase = sprintf('voxelweights_p_%s%s', threshnames{i}, revstr);
            
            montagebase = sprintf('thresh_p_%s_%3.0fvox%s', threshnames{i}, sigvox(i), revstr);
            
            disp(outbase)
            
            % Save weight image
            % ----------------------------------------------
            outname = [outbase '.img'];
            wthr = w;
            wthr(p > thresh(i)) = 0;
            iimg_reconstruct_vols(wthr, volInfo, 'outname', outname);
            
            % ***Note: this can fail if no variability in Z values, i.e. if
            % few booot samples are used. should fix...***
            
            % NEW MONTAGE blobs
            % -------------------------------
            obj = removeblobs(obj);
            
            zthr = sign(w) .* norminv(1 - p');
            zthr(p > thresh(i)) = 0;
%             wpos = wthr; wpos(wpos < 0) = 0;
%             wneg = wthr; wneg(wneg > 0) = 0;
            zpos = zthr; zpos(zpos < 0) = 0;
            zneg = zthr; zneg(zneg > 0) = 0;
            clpos = iimg_indx2clusters(zpos, volInfo);
            clneg = iimg_indx2clusters(zneg, volInfo);
            clpos(cat(1, clpos.numVox) < 3) = [];
            clneg(cat(1, clneg.numVox) < 3) = [];

            obj = addblobs(obj, clpos, 'maxcolor', [1 1 0], 'mincolor', [1 .3 0]);
            obj = addblobs(obj, clneg, 'maxcolor', [0 0 1], 'mincolor', [.7 0 1]);
            
            mb = montagebase; mb(mb == '_') = ' ';
            axes(obj.montage{1}.axis_handles(2))
            title([mb ' k=3 contig'], 'FontSize', 18);

            obj = legend(obj);
            
            scn_export_papersetup(500); saveas(gcf, [montagebase '_new_axial.png']);
            
            diary bootstrap_weightmap_output.txt
            fprintf('\n===============================\n')
            
            disp([mb 'k=3 contig'])
            fprintf('===============================\n')
            fprintf('POSITIVE EFFECTS\n_______________________________\n')
            
            cluster_table(clpos, 0 , 0);
            
            fprintf('\nNEGATIVE EFFECTS\n_______________________________\n')
            
            cluster_table(clneg, 0 , 0);
% %             
% %             % Orthviews
% %             % ----------------------------------------------
% %             cl = iimg_indx2clusters(wthr, volInfo);
% %             cluster_orthviews(cl)
% %             spm_orthviews_hotcool_colormap(w, [prctile(w, 50)]);
% %             cm = get(gcf, 'Colormap');
% %             
% %             % adjust color map
% %             %isgray = all(cm - repmat(mean(cm, 2), 1, 3)  < .01, 2)
% %             cm2 = [cm(1:64, :); negcm(end:-2:1, :); poscm(1:2:end, :)];
% %             colormap(cm2)
% %             
% %             % Montages from orthviews
% %             % ----------------------------------------------
% %             overlay = which('SPM8_colin27T1_seg.img');
% %             
% %             slices_fig_h = cluster_orthviews_montage(6, 'axial', overlay, 0, 'onerow');
% %             %colormap(cm); h = findobj(gcf, 'Type', 'Text'); delete(h)
% %             scn_export_papersetup(500); saveas(gcf, [montagebase '_axial.png']);
% %             
% %             slices_fig_h = cluster_orthviews_montage(6, 'sagittal', overlay, 0, 'onerow');
% %             %colormap(cm); h = findobj(gcf, 'Type', 'Text'); delete(h)
% %             scn_export_papersetup(500); saveas(gcf, [montagebase '_sagittal.png']);
% %             
% %             slices_fig_h = cluster_orthviews_montage(6, 'coronal', overlay, 0, 'onerow');
% %             %colormap(cm); h = findobj(gcf, 'Type', 'Text'); delete(h)
% %             scn_export_papersetup(500); saveas(gcf, [montagebase '_coronal.png']);
% %             
            
        end
        
        % -----------------------------------------------
        % Figure - surface, varying threshold
        % -----------------------------------------------
        
        
        for i = 1:length(thresh)
            
            if ~sigvox(i), continue, end  % skip if no vox
            
            outname = sprintf('thresh_p_%s_%3.0fvox%s', threshnames{i}, sigvox(i), revstr);
            
            disp(outname)
            
            %outbase = 'voxelweights_p002';
            %outname = [outbase '.img'];
            
            wthr = w; wthr(p > thresh(i)+eps) = 0;
            cl = iimg_indx2clusters(wthr, volInfo);
            
            create_figure('Brain_Surface', 2, 2);
            
            
            han1 = addbrain('hires');
            set(han1, 'FaceAlpha', 1);
            cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han1);
            
            %replace with 'hires' and re-run to do whole cortical surface
            create_figure('Brain_Surface', 2, 2, 1);
            subplot(2, 2, 2);
            han = addbrain('hires right'); set(han, 'FaceAlpha', 1);
            
            drawnow
            cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han);
            set(han, 'FaceAlpha', 1);
            axis image; lightRestoreSingle;
            lighting gouraud
            
            create_figure('Brain_Surface', 2, 2, 1);
            subplot(2, 2, 3);
            han = addbrain('hires left'); set(han, 'FaceAlpha', 1);
            axis image; lightRestoreSingle;
            lighting gouraud
            drawnow
            cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han);
            
            create_figure('Brain_Surface', 2, 2, 1);
            subplot(2, 2, 4);
            han = addbrain('limbic'); set(han, 'FaceAlpha', 1);
            axis image; lightRestoreSingle;
            lighting gouraud
            drawnow
            cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han);
            
            set(han(end), 'FaceAlpha', .2)
            han2 = addbrain('brainstem');
            cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han);
            view(135, 15)
            
            create_figure('Brain_Surface', 2, 2, 1);
            subplot(2, 2, 1)
            han = findobj(gca,'Type','patch'); set(han, 'FaceAlpha',1)
            axis image
            
            scn_export_papersetup(600); saveas(gcf, [outname '_Surface1'], 'png');
            subplot(2, 2, 1);
            view(135, 20); lightRestoreSingle;
            saveas(gcf, [outname '_Surface2'], 'png');
            view(225, 20); lightRestoreSingle;
            saveas(gcf, [outname '_Surface3'], 'png');
            view(90, 3); lightRestoreSingle;
            subplot(2, 2, 4);
            view(135, 0); lightRestoreSingle;
            saveas(gcf, [outname '_Surface4'], 'png');
            subplot(2, 2, 1);
            view(270, 3); lightRestoreSingle;
            subplot(2, 2, 4);
            view(235, 0); lightRestoreSingle;
            saveas(gcf, [outname '_Surface5'], 'png');
            subplot(2, 2, 1);
            view(180, -90); lightRestoreSingle;
            saveas(gcf, [outname '_Surface6'], 'png');
            
        end % threshold loop


    end

    % =======================================================================
    % =======================================================================
    % =======================================================================


end % main function



% =======================================================================
% =======================================================================
% =======================================================================
% =======================================================================

% --------------------------------
% Sub-functions
%
% --------------------------------

% =======================================================================
% =======================================================================
% =======================================================================
% =======================================================================

function [v, train_dat, test_dat] = do_pcsquash(train_y, train_dat, test_dat, num_dims, doplssquash)
    % PCA or PLS squash, returns train_dat (scores) and v (weight vectors)
    % and test_dat, with applied weights
    % test_dat is used only to multiply by v to prepare for testing

    % may have to adjust num_dims depending on size of holdout set
    num_dims = min(num_dims, size(train_dat, 1) - 1);

    if doplssquash
        %[v, train_dat] = plssquash(train_dat, train_y, 'num_dims', num_dims, 'noplot');

        [T,P,W,Wstar,U,b,C,Bpls, v, Xhat,Yhat,R2X,R2Y] = PLS_nipals(train_dat,train_y, num_dims);

        % was returning singles sometimes...
        v = double(v);

        % re-do train_dat by taking PLS weighting, [ones train_dat] * v (Bpls_star) = fit
        train_dat = Yhat;
        test_dat = [1 test_dat] * v;

        %     create_figure('tmp1'); plot(Yhat, train_y, 'ko'); drawnow
        %     xlabel('Yhat from PLS'); ylabel('Actual yhat');

    else
        [v, scores] = princomp(train_dat, 'econ');  % scores = train_dat, up to scaling factor!! can't use scores from princomp-diff scaling
        train_dat = train_dat * v;
        v = v(:, 1:num_dims);  % eigenvectors
        train_dat = train_dat(:, 1:num_dims); % train_dat now becomes the scores

        test_dat = test_dat * v; %(:, 1:num_dims(s));  % get scores for missing test subj
    end

end

% ================================
% ================================

function [betas, vox_weights, varargout] = do_fit(fit_method, train_y, train_dat, pcsquash, v, nvox, regparams)
    % subjbetas{s} is a list of voxel weights for one subject, typically
    % (in a linear model) (nvox+1) x 1, where +1 refers to the intercept parameter.
    % in the case of functional mediation, this could be nvox x tpoints.
    %
    % fits are the predicted outcome (y) values, given the model parameter estimates
    % and known information for the left-out observation.  In the functional mediation case,
    % a different method for producing fits is necessary.
    out = [];
    if nargin < 7, regparams = []; end

    switch fit_method
        case 'ols'
            Xs = [ones(size(train_dat, 1), 1) train_dat];
            betas = pinv(Xs) * train_y;

        case 'ridge'
            if isempty(regparams)
                shrinkage_param = 0; % ols
            else
                shrinkage_param = regparams(1);
            end
            Xs = train_dat; % with scaling off, constant is automatically added
            betas = ridge(train_y, Xs, shrinkage_param, 0);

            if pcsquash
                vox_weights = v * betas;
            else
                vox_weights = betas;
            end

        case 'robust'
            Xs = train_dat;
            betas = robustfit(Xs, train_y);

            if pcsquash
                vox_weights = v * betas;
            else
                vox_weights = betas;
            end
            
        case 'bestsubsets'

            if ~pcsquash
                error('Best subsets will not work with over-identified model, so use pcsquash');
            end

            Xs = train_dat;
            if size(Xs, 2) >= size(Xs, 1)-1
                warning('AIC best subsets will not work well with full redundancy.');
            end

            %Xs = Xs(:, 1:end-3); % remove last PC to prevent over-identification (redundancy). AIC will select all predictors with redundant model.
            [wh_predictors, betas] = regress_best_subsets_ga(Xs, train_y);

             
        case 'lasso'
            % Note: train_dat should not have intercept; included automatically
            wh_trace = regparams;
            if size(train_dat, 2) == 1, error('LASSO will not work right with only one predictor variable'); end

            out = lasso(train_y, train_dat);

            %out = lasso_selection(out, 'bic'); % requires mods to original
            %functions; doesn't work for me

            if isempty(out.beta)
                % this is the 'intercept-only' model, where we have
                % information only about the mean (cross-validated)
                % will not work with lasso as implemented though...
                wh_trace = [];
                betas = [];  %out.intercept(wh_trace);

            else
                if isempty(regparams) % if empty, choose OLS
                    wh_trace = size(out.beta, 1);
                else
                    % we have a parameter entered - find which element
                    % corresponds to the input lambda value
                    diffvec = abs(regparams - out.lambda);
                    wh_trace = find(diffvec == min(diffvec));
                    wh_trace = wh_trace(1);

                end
                betas = [out.intercept(wh_trace) out.beta(wh_trace, :)]';
            end

        case 'logistic'

            betas = glmfit(train_dat, [train_y ones(size(train_y))], 'binomial', 'link', 'logit');

        case 'logistictrain'
            [models] = classifierLogisticRegression( train_dat, train_y, [] );
            betas = models{3}(:, 1);  % intercept seems to be added as first predictor

        otherwise error('Unknown fit method')
    end

    if nargout > 2, varargout{1} = out; end

    % Vox weights
    % ------------
    switch fit_method
        case {'ols', 'ridge', 'robust', 'bestsubsets', 'logistic', 'logistictrain'}
            if pcsquash
                vox_weights = v(1:nvox, :) * betas(2:end);
            else
                vox_weights = betas(2:end);
            end

        case 'lasso'

            if pcsquash
                vox_weights = v(1:nvox, :) * out.beta(wh_trace, :)';
            else
                vox_weights = out.beta(wh_trace, :)';
            end

        otherwise error('Unknown fit method')
    end


end % function

% ================================
% ================================



% --------------------------------------------------
% Inner cross-val: Choose best ridge param
% --------------------------------------------------

function  best_paramval = inner_xval_optimize(fit_method, train_y, train_dat, pcsquash, doplot, num_dims, dochoose_ndims, doplssquash, verbose)

    if verbose
        fprintf('\nOptimizing params using inner x-val');
    end

    switch fit_method
        case {'ridge'}
            bounds = [0, 10*size(train_dat, 2)];

        case {'lasso'}
            bounds = [0, .5]; % *****???what should it be?

        otherwise, error('xval_reg:NotImplemented',sprintf('Fit method ''%s'' not compatible with dochoose_regparams', fit_method));
    end

    if pcsquash, pcastr = 'pcsquash';  else pcastr = 'noverbose'; end  % set PCA string
    if pcsquash && doplssquash, plsstr = 'pls';  else plsstr = 'noverbose'; end  % set PLS option (not tested)

    t0 = clock;

    % paramval could be a ridge or LASSO shrinkage value
    fit_fcn = @(paramval) xval_given_param(paramval, fit_method, train_y, train_dat, pcastr, plsstr, num_dims);

    best_paramval = fminbnd(fit_fcn, bounds(1), bounds(2));

    if verbose
        fprintf(': Done: %3.2f sec. ', etime(clock, t0));
        fprintf('Best regularization parameter value = %3.2f\n', best_paramval)
    end

end  % inner_xval_optimize_ridge

% PE-generating function for ridge

function pe = xval_given_param(paramval, fit_method, Y, data, pcastr, plsstr, num_dims)
    % paramval could be a ridge or LASSO shrinkage value
    STATS = xval_regression_multisubject(fit_method, {Y}, {data}, 'regparams', paramval, 'noverbose', 'holdout_method', 'balanced4', pcastr, plsstr, 'num_dims', num_dims);
    pe = STATS.pred_err;
end
