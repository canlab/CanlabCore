function STATS = xval_regression_multisubject(fit_method, Y, X, varargin)
    % STATS = xval_regression_multisubject(fit_method, Y, X, varargin)
    %
    % CROSS-VALIDATED Regression
    % Leave-one observation out, predict outcomes for each missing holdout_set.
    % OR define other types of holdout sets, as described below
    %
    % Y = outcome, N x 1, enter for each separate dataset (e.g., subject) in cells {}
    % X = predictors, N x variables, enter for each separate dataset (e.g., subject) in cells {}
    % fit_method = 'lasso' 'ols' 'ridge' or 'robust'
    %
    % optional inputs:
    % 'pca', PCA data reduction first
    % 'ndims', dims to save from PCA
    % 'variable', retain max dims for each subject
    % case {'cov', 'covs'}, cov = varargin{i+1};  NOTE:covariates are added to training data if entered.
    %
    % {'lassopath', 'ridgek', 'regparams'}, regparams = varargin{i+1};
    %  - For LASSO: This is the lambda reg param, should be between 0 and 1
    %
    % {'noverb', 'noverbose'}, verbose = 0;
    % {'lowverb', 'loverb', 'lowverbose'}, verboseL = 1; verbose = 0;
    % 'nested_choose_ndims', dochoose_ndims = 1;
    %
    % 'optimize_regularization', optimize LASSO lambda or ridge parameter
    % using inner x-val loop with balanced4 selection and least squares
    % nonlinear fit
    %
    % 'holdout_method',     followed by the name of a holdout method:
    % 'loo'                     leave-one-out
    % 'categorical_covs'    balanced on categorical covariates
    % 'l4o_covbalanced'     leave-4-out, balanced on combo of Y and continuous
    %                           covariates (logistic regression/propensity score based method)
    % 'balanced4'           4-fold, balanced on Y (every 4th element of sorted Y)
    % Or CUSTOM holdout set:Enter integer list of which obs belong to which holdout set. 
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
    % Tor Wager, 3/17/09

    % Programmers' notes
    % June 19, 2011: minor change to LASSO lambda selection when entering
    % regparam
    %
    % July 31, 2011: change to weight map: use whole sample for
    % weights/betas
    %
    % Feb 6, 2012: Minor bug fix in multisubject mode
    %
    % Nov 10, 2012: Tor Wager: Updated name of lasso to lasso_rocha.m for
    % matlab compatibilty.  Forced double-format for vars to avoid
    % single-format problems.
    
    % Set defaults
    % -----------------------------------------------------------
    nested_setdefaults();
    docenterrowsX; dosave;

    % Variable dimensions: retain max possible for each subject
    % If choose ndims is selected, then pick the number here
    % Based on size of input data only (not optimized; full dims)
    % -----------------------------------------------------------
    nested_choose_ndims();

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
        holdout_set = nested_select_holdout_set();
        STATS.INPUTS.holdout_set{s} = holdout_set;

        checkrowdependence(X{s}, verbose || verboseL);  
        % rows are assumed to be independent; if they are not, can get
        % perfect cross-validated predictions of even random outcomes
        
        % --------------------------------------------------------------------------------
        % Cross-validation folds
        % --------------------------------------------------------------------------------

        % THIS COULD BE  PARFOR
        for wh_fold = 1:length(holdout_set)  %holdout_set = 1:length(Y{s})

            if verbose || verboseL, fprintf('\b\b\b%3.0f', wh_fold); end

            % Parse into training and test data for this fold
            % Defines: train_y (outcome) and train_dat (predictors),
            %              test_dat.  test outcome not saved.
            % covariates are added to training data if entered.
            nested_select_training_test_data();
            
            % Parameter optimization by inner cross-validation
            % Inner cross-validation could take a long time
            % --------------------------------------------------------------------------------

            if dochoose_ndims
                % update num_dims(s), update regparams***incorporate this
                % into single inner-xval loop
                [LASSO, regparams, best_ndims] = inner_xval_lasso(fit_method, train_y, train_dat, pcsquash, doplot, num_dims, dochoose_ndims, doplssquash);

            end

            if dochoose_regparams
                % Replace regparams with optimized regparams
                regparams = inner_xval_optimize(fit_method, train_y, train_dat, pcsquash, doplot, num_dims, dochoose_ndims, doplssquash, verbose);

                STATS.all_reg_hyperparams(wh_fold) = regparams;
            end

            % Dim reduction
            % --------------------------------------------------------------------------------
            if pcsquash
                [v, train_dat, test_dat] = do_pcsquash(train_y, train_dat, test_dat, num_dims(s), doplssquash);
            end

            % Fit
            % --------------------------------------------------------------------------------
            % subjbetas{s} is a list of voxel weights for one subject, typically
            nvox = size(X{s}, 2);  % original voxels, not including covs or NaN voxels (will add in NaNs later)
            
            % This is for the version where covs are always added after
            % pcsquash; but if covs are not meaningful, can increase
            % error!
%             ncovs = size(train_covs, 2); % needed to get voxel weights only
%             [subjbetas{s}, STATS.vox_weights(:, wh_fold)] = do_fit(fit_method, train_y, [train_covs train_dat], pcsquash, v, nvox, regparams, ncovs);

            [subjbetas{s}, STATS.vox_weights(:, wh_fold)] = do_fit(fit_method, train_y, train_dat, pcsquash, v, nvox, regparams);

            switch fit_method
                case {'logistic', 'logistictrain'}
                    eta  = [1 test_dat] * subjbetas{s};
                    fit(holdout_set{wh_fold}, 1) = 1 ./ (1 + exp(-eta));

                otherwise
                    % pred for the left-out observation
                    if length(subjbetas{s}) == 1
                        % this is the 'intercept-only' model, where we have
                        % information only about the mean (cross-validated)
                        fit(holdout_set{wh_fold}, 1)  = 1 * subjbetas{s};
                    else
                        %fit(holdout_set{wh_fold}, 1)  = [ones(size(test_dat, 1), 1) test_covs test_dat] * subjbetas{s};
                        fit(holdout_set{wh_fold}, 1)  = [ones(size(test_dat, 1), 1) test_dat] * subjbetas{s};
                    end
            end

        % plot apparent fit
%         if verbose
%             if wh_fold ==1, create_figure('Intercepts'); end
%             create_figure('Training data', 1, 2); plot([ones(size(train_dat, 1), 1) train_dat] * subjbetas{s}, train_y, 'ko', 'MarkerFaceColor', [.5 .5 .5]);
%             hold on; plot(fit(holdout_set{wh_fold}, 1), Y{s}(holdout_set{wh_fold}), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
%             create_figure('Intercepts', 1, 1, 1); hold on; plot(holdout_set{wh_fold}, subjbetas{s}(1), 'ks', 'MarkerFaceColor', 'k'); title('Intercept'); drawnow
%         end
        
% Diagnostic: Is apparent fit related to covariates?
% train_covs = cov_val{1}; train_covs(holdout_set{wh_fold}, :) = [];
% [b, dev, statsss] = glmfit(train_covs, [ones(size(train_dat, 1), 1) train_dat] * subjbetas{s});
% glm_table(statsss);

        my_intercepts{s}(wh_fold, 1) = subjbetas{s}(1);
        
        end % xval

        if verbose || verboseL, toc, end

        % ---------------------------------------------------------------------
        % end cross-val loop
        % ---------------------------------------------------------------------

        subjfit{s} = fit;

        nested_output_metrics_and_plots();
        
        % ---------------------------------------------------------------------
        % Get weights and betas across ALL observations
        % Tor added 7/31/11
        % ---------------------------------------------------------------------
        
        if verbose
            disp('Getting final variable weights across all observations');
        end
        
        weight_y = Y{s};
        
        if ~isempty(cov_val)
            weight_dat = [X{s} cov_val{s}];
        else
            weight_dat = X{s};
        end
        
        % Dim reduction
        % --------------------------------------------------------------------------------
        if pcsquash
            [v, weight_dat] = do_pcsquash(weight_y, weight_dat, weight_dat, num_dims(s), doplssquash);
        end
        
        % subjbetas are final weights after OLS
        [subjbetas{s}, STATS.mean_vox_weights(:, s)] = do_fit(fit_method, weight_y, weight_dat, pcsquash, v, nvox, regparams);
        
        if pcsquash
            % add NaNs back in to preserve voxel order
            STATS.mean_vox_weights(:, s) = naninsert(nanvox{s}, STATS.mean_vox_weights(:, s));
            
        end

        
    
    end % dataset/subject loop

 
                
    STATS.Y_orig = Y_orig;
    if ~isempty(cov_val)
        STATS.note = 'covariates were added to predictive model'; %, if covs %removed in Y_orig stored here, if covs were entered';
    else
        STATS.note = 'no covariates entered';
    end

    STATS.my_intercepts = my_intercepts;
    
    STATS.devs_from_mean_only_model = devs_from_mean_only_model;
    STATS.devs_from_full_model = devs_from_full_model;
    STATS.var_null = var_null;
    STATS.var_full = var_full;

    % Null model leave-one-out prediction error
    % By predicting based on the mean of OTHER subjects, we've increased
    % the variance
    STATS.pred_err_null = var_null .^ .5;

    STATS.pred_err = pred_err;
    STATS.pred_err_descrip = 'Apparent error/loss: Root mean squared deviation from observed outcomes';
    STATS.var_reduction = rsq;
    STATS.var_reduction_descrip = 'Percent reduction in outcome variance due to model; negative values indicate added variance.';

    STATS.subjfit = subjfit;
    STATS.subjbetas = subjbetas;
    STATS.r_each_subject = rr;
    STATS.r_squared = STATS.r_each_subject .^ 2;
    STATS.r_each_subject_note = 'r value: correlation between predicted and outcome values';
    STATS.r_each_subject_note2 = 'may not be very interpretable, because null model r = -1.0';

    if verbose
        disp(['Mean correlation is: ' num2str(mean(rr))])
        disp(['Mean proportion of original variance explained is: ' num2str(mean(rsq))])
        disp('The number above is based on reduction of original variance;');
        disp('it can be negative if predictors are not helpful because they add noise, increasing the overall variance!')
        disp(' ')
        fprintf('Null model pred. error is %3.2f, and full model is %3.2f\n', STATS.pred_err_null(1), STATS.pred_err(1));
        disp(' ')
    end



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
        docenterrowsX = 0;
        dosave = 1;
        
        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}
                    % Data
                    case {'centerX', 'centerx', 'docenterrowsX'}, docenterrowsX = 1;
                        
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

        if strcmp(fit_method, 'lasso')
            if ~exist('lasso_rocha.m', 'file')
                addpath(lassopath)
            end
        end

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

        % force double 
        X{s} = double(X{s});
        Y_orig{s} = double(Y_orig{s});
        
        if docenterrowsX 
            disp('Row-centering requested: Centering each row of predictors');
            X{s} = scale(X{s}', 1)'; 
        end
        
        fit = NaN * zeros(size(Y{s}));
        tic

        % ---------------------------------------------------------------------
        % initialize optional things
        v = [];

        if verbose || verboseL, fprintf('Fold %3.0f', 0); end

    end

     % =======================================================================   
    function  checkrowdependence(X, verbose)  
        % rows are assumed to be independent; if they are not, can get
        % perfect cross-validated predictions of even random outcomes

        [N, k] = size(X);
        r = rank(X');
        if r < N && verbose, disp('WARNING!!! ROWS ARE DEPENDENT. CROSS-VALIDATION MAY BE INVALID.'); end
%         
%         cc = corrcoef(X');
%         cc = cc - eye(size(cc));
%         [maxcorr, wh] = max(cc(:));
% 
%         vifs = getvif(dat', 1)
%         vifs = getvif(dat', 1); [maxvif, whv] = max(vifs)

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


    function holdout_set = nested_select_holdout_set
        % purpose: return holdout_set variable
        
        nobs_s = length(Y{s});
        
        if ~ischar(holdout_method)
            % assume holdout method entered is actual holdout set, integer
            % list of which obs belong to which test set
            u = unique(holdout_method);
            
            holdout_set = cell(1, length(u));
            for i = 1:length(u)
                holdout_set{i} = find(holdout_method == i);
            end
            
            return
        end
        
        switch lower(holdout_method)
            case 'loo'
                holdout_set = cell(1, nobs_s);
                for i = 1:nobs_s, holdout_set{i} = i; end

            case 'l4o_covbalanced'
                disp('Selecting holdout sets: Balancing on covariates and outcome, and also trying to ensure that each obs is selected equally often.');
                holdout_proportion = 4 ./ nobs_s; % prop for leave-4-out
                nfolds = 40;
                if isempty(cov_val)
                    wh_holdout = xval_select_holdout_set(Y{s}, [], nfolds, holdout_proportion, verbose);
                else
                    wh_holdout = xval_select_holdout_set(Y{s}, cov_val{s}, nfolds, holdout_proportion, verbose);
                end
                holdout_set = cell(1, nfolds);
                for k = 1:nfolds
                    holdout_set{k} = wh_holdout(:, k);
                end

            case 'categorical_covs'
                
                holdout_set = xval_select_holdout_set_categoricalcovs(cov_val{s});
                
            case 'balanced4'
                nfolds = 4;
                holdout_set = cell(1, nfolds);
                [ys, indx] = sort(Y{s});
                for k = 1:nfolds

                    holdout_set{k} = indx(k:nfolds:end);
                    if isempty(holdout_set{k}), error('Holdout set construction error: Dataset too small?'); end

                end

            otherwise error('Unknown holdout method.  See help.');
        end
    end


    % =======================================================================
    
    function nested_select_training_test_data
        % Select training/test data
        % -----------------------------

        train_y = Y{s};
        train_y(holdout_set{wh_fold}) = [];

        if ~isempty(cov_val)
            train_dat = [X{s} cov_val{s}];
        else
            train_dat = X{s};
        end

                    % This is for the version where covs are always added after
            % pcsquash; but if covs are not meaningful, can increase
            % error!
% %         if ~isempty(cov_val)
% %             train_covs = cov_val{s};
% %             test_covs = train_covs(holdout_set{wh_fold}, :);
% %             train_covs(holdout_set{wh_fold}, :) = [];              % leave out the missing observation(s)
% %         else
% %             train_covs = [];
% %             test_covs = [];
% %         end

        test_dat = train_dat(holdout_set{wh_fold}, :);
        train_dat(holdout_set{wh_fold}, :) = [];              % leave out the missing observation(s)


    end % data selection function

    % =======================================================================


    function nested_output_metrics_and_plots


        pefcn = inline('var(y - f) .^ .5', 'y', 'f');
        
        if verbose || verboseL

            create_figure('fitplot', 2, 2);
            plot(Y_orig{s}, 'o-');
            hold on; plot(fit, 'rx-');
            title('Predicted (red), Actual (blue) over obs');

            subplot(2, 2, 2);
            plot_correlation_samefig(fit, Y_orig{s});
            xlabel('Predicted value'); ylabel('Actual value')

            plot([min(Y_orig{s})  max(Y_orig{s})], [min(Y_orig{s})  max(Y_orig{s})], 'b--', 'LineWidth', 2);
            axis tight

            subplot(2, 2, 3);
            hold on; plot(1:length(my_intercepts{s}), my_intercepts{s}, 'ks', 'MarkerFaceColor', 'k');
            title('Intercept');

            subplot(2, 2, 4);
            rescalevals = [0:.1:2];
            for j = 1:length(rescalevals)
                mype(j) = pefcn(Y_orig{s}, rescalevals(j) .* subjfit{s});
            end
            plot(rescalevals, mype, 'k-');
            plot_vertical_line(0);
            axis auto; axis tight;
            title('Pred err as a function of prediction rescaling');
            drawnow

            if ~isempty(cov_val)
                if dosave, diary('xval_regression_multisubject_output.txt'); end

                disp('Test of whether predicted values are related to covariates:');
                [b, dev, statsss] = glmfit(cov_val{1}, subjfit{s});
                glm_table(statsss);

                if dosave, diary off; end
            end

        end % verbose

        myr = corrcoef(subjfit{s}, Y_orig{s}); rr(s) = myr(1,2);
        pred_err(s, 1) = var(Y_orig{s} - subjfit{s}) .^ .5; % almost the same as (sum((Y_orig{s} - subjfit{s}).^2) ./ (length(Y_orig{s}) - 1) ) .^.5

        % Null model leave-one-out prediction error
        % By predicting based on the mean of OTHER subjects, we've increased
        % the variance; we need to adjust
        try
            
        jstat = jackknife(@mean, Y_orig{s});

        catch
            disp('Weird Matlab 2010a jackknife error.  Debug me?!!');
            jstat = Y_orig{s};
        end
        
        devs_from_mean_only_model = Y_orig{s} - jstat;

        devs_from_full_model = Y_orig{s} - subjfit{s};

        var_null(s) = var(devs_from_mean_only_model);
        var_full(s) = var(devs_from_full_model);

        rsq(s, 1) = 1 - var_full(s) ./ var_null(s);

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

            % This is for the version where covs are always added after
            % pcsquash; but if covs are not meaningful, can increase
            % error!
       %function [betas, vox_weights, varargout] = do_fit(fit_method, train_y, train_dat, pcsquash, v, nvox, regparams, ncovs)
     
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
            %betas = pinv(Xs) * train_y; 
            % this produces DIFFERENT results from standard with v >> N
            % better weight stability, but also apparent positive bias (in
            % weights?)
            
            %betas = inv(Xs'*Xs) * Xs' * train_y;
            betas = pinv(Xs) * train_y;
            
            % changed back, Oct 2011, because inv too slow, and pinv should
            % be more stable overall.
            
        case 'ridge'
            if isempty(regparams)
                shrinkage_param = 0; % ols
            else
                shrinkage_param = regparams(1);
            end
            Xs = train_dat; % with scaling off, constant is automatically added
            betas = ridge(train_y, Xs, shrinkage_param, 0);


        case 'robust'
            Xs = train_dat;
            betas = robustfit(Xs, train_y);


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
            wh_trace = regparams; % empty, or takes input lambda value
            if size(train_dat, 2) == 1, error('LASSO will not work right with only one predictor variable'); end

            if ~exist('lasso_rocha.m', 'file')
                error('THE LASSO USED HERE HAS BEEN RENAMED LASSO_ROCHA FOR MATLAB COMPATIBILITY. LASSO_ROCHA.M MUST BE ON YOUR MATLAB PATH.');
            end
                
            out = lasso_rocha(train_y, train_dat);

            %out = lasso_selection(out, 'bic'); % requires mods to original
            %functions; doesn't work for me
            
            if isempty(out.beta)
                % this is the 'intercept-only' model, where we have
                % information only about the mean (cross-validated)
                % will not work with lasso as implemented though...
                wh_trace = [];
                betas = [];  %out.intercept(wh_trace);

            else
                if isempty(wh_trace) % if empty, choose OLS
                    wh_trace = size(out.beta, 1);
                else
                    % we have a parameter entered - find which element
                    % corresponds to the input lambda value
                    notok = all(out.beta == 0, 2); % can return all invalid values for b
                    diffvec = abs(regparams - out.lambda);
                    diffvec(notok) = Inf;
                    wh_trace = find(diffvec == min(diffvec));
                    wh_trace = wh_trace(end);

%                     create_figure('trace', 2, 1);
%                     plot(out.lambda, out.beta);
%                     h = plot_vertical_line(regparams);
%                     h2 = plot_vertical_line(out.lambda(wh_trace));
%                     set(h2, 'Color', 'r');
%                     subplot(2, 1, 2);
%                     plot(out.beta(wh_trace, :)); title('Betas')
%                     drawnow
%                     pause(1);
                end
                
                    % re-fit with OLS, using only non-zero elements
                    % see Hastie et al., "Elements" book, p. 92
                    
                    wh_in_model = logical([1 (out.beta(wh_trace, :) ~= 0)]');
                    
                    Xs = [ones(size(train_dat, 1), 1) train_dat(:, wh_in_model(2:end))];

                    betatmp = pinv(Xs) * train_y;
                    
                    betas = zeros(size(wh_in_model));
                    betas(wh_in_model) = betatmp;
                    
                %betas = [out.intercept(wh_trace) out.beta(wh_trace, :)]';
            end

        case 'logistic'

            betas = glmfit(train_dat, [train_y ones(size(train_y))], 'binomial', 'link', 'logit');

            % plot apparent
            % %                 eta  = [ones(size(train_dat, 1), 1) train_dat] * betas;
            % %                 fit = 1 ./ (1 + exp(-eta));
            % %                 create_figure('tmp');
            % %                 plot(train_dat, train_y, 'ko');
            % %                 vals = [-1.5:.1:1.5]';
            % %                 etafit = [ones(length(vals), 1) vals] * betas;
            % %                 fitcurve = 1 ./ (1 + exp(-etafit));
            % %                 plot([-1.5:.1:1.5], fitcurve, 'k');
            % %                 set(gca, 'XLim', [min(train_dat) max(train_dat)]);
            % %                 etazero = -betas(1) ./ betas(2);  % point at which p(lie) in logistic model = .5
            % %                 han = plot_vertical_line(etazero);
            % %                 set(han, 'LineStyle', ':');
            % %                 drawnow

            % % %     case 'logisticpls'
            % % %         now done in plssquash...
            % % %          [T,P,W,Wstar,U,b,C,Bpls,Bpls_star,Xhat,Yhat,R2X,R2Y] = PLS_nipals(train_dat,train_y, 10);
            % % %          betas was Bpls_star; [1 test_x] * Bpls_star = fit = trial score,
            % % %          which is X in the logistic model below:
            % % %
            % % %         Xtrain = [ones(size(train_dat, 1), 1) train_dat] * Bpls_star;  % this gives us the PLS-optimized features
            % % %         betas = glmfit(Xtrain, [train_y ones(size(train_y))], 'binomial', 'link', 'logit');
            % % %
            % % %         fit_eta = betas(1) + betas(2) * ([1 test_data]*Bpls_star)
            % % %         fit = 1 ./ (1 + exp(-fit_eta));

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
                %vox_weights = v(1:nvox, :) * betas(2+ncovs:end);
                vox_weights = v(1:nvox, :) * betas(2:end);
            else
                vox_weights = betas(2:end);
            end

        case 'lasso'

            % This is for the version where covs are always added after
            % pcsquash; but if covs are not meaningful, can increase
            % error!
            %             if pcsquash
            %                 vox_weights = v(1:nvox, :) * out.beta(wh_trace, ncovs+1:end)';
            %             else
            %                 vox_weights = out.beta(wh_trace, ncovs+1:end)';
            %             end

            % NOTE: tor edited to use re-fit OLS solution in case of
            % penalization.
            
            if pcsquash
                vox_weights = v(1:nvox, :) * betas(2:end); % out.beta(wh_trace, :)';
            else
                vox_weights = betas(2:end); % out.beta(wh_trace, :)';
            end

        otherwise error('Unknown fit method')
    end


end % function

% ================================
% ================================

% % % function [LASSO, wh_trace, best_ndims] = inner_xval_lasso(fit_method, train_y, train_dat, pcsquash, doplot, num_dims, dochoose_ndims, doplssquash)
% % %
% % % % Inner cross-validation; this could take a long time
% % % % -----------------------------
% % %
% % % % update num_dims(s), update wh_trace
% % %
% % % n_inner_obs = length(train_y);
% % % nvox = size(train_dat, 2);  % original voxels, not including covs or NaN voxels (will add in NaNs later)
% % %
% % % best_ndims = size(train_dat, 2);  % initialize num_dims value
% % %
% % % % choose dims to test -- either variable, or single value
% % % % --------------------------------------------------------------------
% % % if pcsquash && dochoose_ndims
% % %     dims_to_test = 2:num_dims; % placeholder; select all features
% % % else
% % %     % run once, select all features
% % %     dims_to_test = num_dims;
% % % end
% % %
% % % % choose lasso shrinkage parameter with inner cross-validation loop
% % % % -------------------------------------------------------------------
% % % if strcmp(fit_method, 'lasso')
% % %     if pcsquash
% % %         % 2nd train_dat(1,:)is test dat, but irrelevant here-- we just need dims
% % %         [v, X2, tmp_train_dat] = do_pcsquash(train_y, train_dat, train_dat(1,:), dims_to_test(end), doplssquash);
% % %         clear v tmp_train_dat
% % %     end
% % %     out = lasso(train_y, X2); % to get penalty values
% % %     penalty_values = out.penalty;  %linspace(0, max(out.penalty), 10);  % resolution: 10
% % %
% % %     fit = zeros(n_inner_obs, length(dims_to_test), length(penalty_values));
% % %
% % % end
% % %
% % % % Inner cross-validation
% % % % --------------------------
% % % for ii = 1:n_inner_obs
% % %
% % %     % fprintf('\b\b\b%3.0f', holdout_set);
% % %
% % %     % Select data
% % %     Y2 = train_y;
% % %     Y2(ii) = [];
% % %     X2 = train_dat;
% % %     X2(ii, :) = [];
% % %     test_dat = train_dat(ii, :);
% % %
% % %     % Dim reduction: with max number
% % %     % We will use selected numbers of these
% % %     % Must re-do PCA in inner loop to avoid bias
% % %     % -----------------------------
% % %     if pcsquash
% % %         [v, X2, test_dat] = do_pcsquash(Y2, X2, test_dat, dims_to_test(end), doplssquash);
% % %     end
% % %
% % %     for jj = 1:length(dims_to_test)
% % %         % Dimension selection loop
% % %
% % %         % Fit
% % %         % -----------------------------
% % %         % we may also have a lasso shrinkage (penalty) parameter to choose
% % %         % this gives us the whole lasso trace, though
% % %         wh_trace = [];
% % %         [subjbetas{ii,jj}, tmp_vox_weights, out2] = do_fit(fit_method, Y2, ...
% % %             X2(:, 1:dims_to_test(jj)), pcsquash, v(:, 1:dims_to_test(jj)), ...
% % %             nvox, wh_trace);
% % %
% % %         clear tmp_vox_weights
% % %
% % %         % pred for the left-out observation
% % %
% % %
% % %         if strcmp(fit_method, 'lasso')
% % %             % get trace from out struct
% % %             % penalty values vary across validation loop iterations
% % %             % fit = fit values for each possible choice of penalty
% % %             % test_dat is component scores if pcsquash, otherwise dims to
% % %             % test is all
% % %             fit(ii, jj, 1) = out2.intercept(end)' + test_dat(1:dims_to_test(jj)) * out2.beta(end, :)';
% % %
% % %             % lfit: lasso fit for each trace
% % %             lfit = out2.intercept' + test_dat(1:dims_to_test(jj)) * out2.beta';
% % %             % interpolate to standard trace values
% % %             lfit = interp1(out2.penalty, lfit, penalty_values);
% % %
% % %             fit(ii, jj, 1:length(penalty_values)) = lfit;
% % %             % this should go on 3rd dim of fit
% % %
% % %             % trials x penalty choices
% % %             %             fit_inner_xval{jj}(ii, :) = interp1(out2.penalty, fitiijj, penalty_values);
% % %             %             err_inner_xval{jj}(ii, :) = train_y(ii) - fit_inner_xval{jj}(ii, :);
% % %         else
% % %
% % %             % this is for the OLS solution (max trace value)
% % %             if length(subjbetas{ii,jj}) == 1
% % %                 % this is the 'intercept-only' model
% % %                 fit(ii, jj, 1)  = 1 * subjbetas{ii,jj};
% % %             else
% % %                 fit(ii, jj, 1)  = [1 test_dat] * subjbetas{ii,jj};
% % %             end
% % %         end
% % %
% % %
% % %
% % %     end % Dims loop
% % % end % ii inner xval loop
% % %
% % % for j = 1:length(dims_to_test)
% % %
% % %     for k = 1:size(fit, 3)
% % %         vec = squeeze(fit(:, j, k));
% % %         cc = corrcoef(vec, train_y);
% % %         r(j, k) = cc(1,2);
% % %         pe(j, k) = sqrt(sum((train_y - vec) .^2));
% % %     end
% % %
% % % end
% % %
% % % [wh_ndims, wh_penalty] = find(pe == min(pe(:)));
% % % wh_ndims = wh_ndims(1);
% % % wh_penalty = wh_penalty(1);
% % %
% % % best_penalty = penalty_values(wh_penalty);
% % % best_dims = dims_to_test(wh_ndims);
% % % best_pe = pe(wh_ndims, wh_penalty);
% % %
% % % if doplot
% % %     create_figure('Error by dims and penalty values');
% % %
% % %     [Xf, Yf] = meshgrid(penalty_values, dims_to_test);
% % %     surf(Xf, Yf, pe)
% % %     xlabel('Penalty values')
% % %     ylabel('Dimensions')
% % %     view(135, 30);
% % %     hold on;
% % %     plot3(best_penalty, best_dims, best_pe, 'ko','MarkerFaceColor', 'r', 'MarkerSize', 8);
% % %
% % % end
% % %
% % % LASSO = struct('best_penalty', best_penalty, 'best_dims', best_dims, 'best_pe', best_pe, ...
% % %     'wh_ndims', wh_ndims, 'wh_penalty', wh_penalty);
% % %
% % % end % function


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
            bounds = [0, 1]; % bounds on lambda
            
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
