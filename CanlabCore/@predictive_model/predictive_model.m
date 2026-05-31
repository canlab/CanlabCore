classdef predictive_model
    % predictive_model  Object representing a multivariate predictive model.
    %
    % Canonical output of CanlabCore's predictive-modelling entry points
    % (xval_SVM, xval_SVR, xval_discriminant_classifier, xval_cross_classify,
    % xval_regression_multisubject, xval_lasso_brain, xval_ridge_brain,
    % xval_bestsubsets_brain, fmri_data.predict, ...).
    %
    % TOP-LEVEL PROPERTIES (post-consolidation)
    %
    %   HYPERPARAMETERS (set by user/constructor)
    %     algorithm, task, modeloptions, random_state, standardize,
    %     use_parallel, cv, scorer, nboot, nperm, do_calibrate,
    %     Y_name, X_name, class_labels
    %
    %   DATA (top-level, was in the inputs category)
    %     Y                 outcome vector actually fit (after bad-case removal)
    %     id                grouping vector (e.g., subject id)
    %     omitted_cases     logical, length = original Y, true = removed
    %                       by the pre-fit data-quality check (NaN in X or Y,
    %                       or any other "bad case" the wrapper detected).
    %     omitted_features  logical, length = original feature count,
    %                       true = removed (all-NaN columns, zero-variance
    %                       columns, etc.).
    %     inputParameters   struct: pcsquash, num_dims, holdout_method,
    %                       regparams, lambda, Alpha, snapshot of modeloptions,
    %                       etc. Consolidates legacy INPUTS, inputs, inputOptions.
    %
    %   FIT METADATA
    %     fit_type          char: 'crossval' | 'insample' | 'test'.
    %                       Describes the provenance of BOTH `yfit` and `weights`.
    %                       Set by the function that populates them; for
    %                       cross-validated wrappers this is 'crossval'.
    %
    %   FITTED STATE (categorised sub-structs)
    %     fitted_values     .yfit  .scores  .score_type ('distance'|'probability')
    %                       .scores_within_id .Y_within_id .scorediff
    %                       .high_vs_low_scores_within_id .subjfit
    %                       .predictions .Y_per_fold
    %     weights           .w (main coefficient vector)
    %                       .w_perfold (matrix nvox x nfolds; was vox_weights)
    %                       .w_bootstrap (struct; was VOXWEIGHTS)
    %                       .mean_vox_weights .my_intercepts .subjbetas .weight_obj
    %                       .boot_w .boot_w_ste .boot_w_mean      (merged from weight_stats)
    %                       .z .p .fdr_thr .fdr_sig .thresh_fdr   (merged from weight_stats;
    %                                                             renamed wZ->z, wP->p,
    %                                                             wP_fdr_thr->fdr_thr,
    %                                                             boot_w_fdrsig->fdr_sig,
    %                                                             w_thresh_fdr->thresh_fdr)
    %     error_metrics     Each entry is a (value, descrip) tuple:
    %                       error_metrics.pred_err = struct('value', 2.4, 'descrip', 'RMSE')
    %                       Top-level Dependent aliases (.pred_err, .crossval_accuracy, ...)
    %                       unwrap .value for arithmetic; the corresponding
    %                       *_descrip aliases unwrap .descrip.
    %     cv_partition      .trIdx .teIdx .nfolds .fold_modeloptions
    %                       .hyperparams_by_fold .cvpartition
    %     ml_model          trained MATLAB model object (verbatim).
    %     fold_models       cell of per-fold trained models.
    %     bootstrap_results full bootstrap output struct.
    %     permutation_results
    %     diagnostics       .ROC_forced_choice .ROC_single_interval
    %                       .mult_obs_within_person .all_reg_hyperparams .phi
    %     cross_classify    cross-classification output (.stats1 .stats2
    %                       .test_results .crosstestfit .cvoutput .test_Y .all)
    %     legacy_extras     catch-all for fields without a categorised home.
    %     note              cell array of strings, appended to (not overwritten)
    %                       as wrappers add commentary. Consolidates legacy
    %                       note + r_each_subject_note + r_each_subject_note2.
    %     accfun            objective function handle (top-level).
    %     history           cell of strings (image_vector convention).
    %     removed_voxels    logical (image_vector convention).
    %
    % LEGACY DEPENDENT ALIASES (read-only, full back-compat)
    %     Every flat field name that user code currently reads still works:
    %     Y_orig, trueLabels, y -> Y
    %     INPUTS, inputs, inputOptions -> inputParameters
    %     vox_weights -> weights.w_perfold
    %     VOXWEIGHTS  -> weights.w_bootstrap
    %     wZ, wP, wP_fdr_thr, boot_w_fdrsig, w_thresh_fdr
    %                 -> weights.z / .p / .fdr_thr / .fdr_sig / .thresh_fdr
    %     r_each_subject_note, r_each_subject_note2 -> indexed entries in .note
    %     pred_err, crossval_accuracy, r_squared, mse, rmse, ... -> .value of the tuple
    %     pred_err_descrip, crossval_accuracy_descrip, ...       -> .descrip of the tuple
    %     classification_d_singleinterval, regression_d_singleinterval, ...
    %                                                            -> d_singleinterval.value
    %     SVMModel, SVRModel, ClassificationModel -> ml_model
    %
    % :Usage:
    % ::
    %     pm = predictive_model();             % empty
    %     pm = predictive_model(S);            % from struct (legacy wrappers)
    %     pm = predictive_model(S, 'noverbose'); % suppress "fields not copied"
    %     pm = predictive_model('algorithm','svm','task','classification', ...
    %                           'modeloptions',{'KernelFunction','linear'});
    %
    % :Examples:
    % ::
    %     S  = xval_SVM(X, Y, id, 'nooptimize', 'norepeats', 'nobootstrap');
    %     pm = predictive_model(S);
    %     pm.fit_type                  % 'crossval'
    %     pm.crossval_accuracy         % numeric (Dependent unwraps .value)
    %     pm.crossval_accuracy_descrip % descrip string (Dependent unwraps .descrip)
    %     pm.weights.z                 % bootstrap z-scores
    %     pm.weights.w_perfold         % per-fold weight matrix
    %
    % :See also:
    %   xval_SVM, xval_SVR, xval_discriminant_classifier, fmri_data.predict
    %
    % -------------------------------------------------------------------------
    %     Author and copyright information:
    %
    %     Copyright (C) 2025  CANlab
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    % -------------------------------------------------------------------------

    properties  % HYPERPARAMETERS
        algorithm       = '';
        task            = '';
        modeloptions    = {};
        random_state    = [];
        standardize     = [];
        use_parallel    = [];
        cv              = [];
        scorer          = [];
        nboot           = [];
        nperm           = [];
        do_calibrate    = [];
        Y_name          = '';
        X_name          = '';
        class_labels    = {};
    end

    properties (SetAccess = protected)  % DATA + FIT METADATA + CATEGORISED FITTED STATE
        Y                    = [];
        id                   = [];
        omitted_cases        = [];
        omitted_features     = [];
        inputParameters      = struct();
        fit_type             = '';

        cv_partition         = struct();
        fitted_values        = struct();
        weights              = struct();
        error_metrics        = struct();
        ml_model             = [];
        fold_models          = {};
        bootstrap_results    = struct();
        permutation_results  = struct();
        diagnostics          = struct();
        cross_classify       = struct();
        legacy_extras        = struct();
        note                 = {};
        accfun               = [];
        history              = {};
        removed_voxels       = [];
    end

    properties (Dependent, SetAccess = private)  % LEGACY READ-ONLY ALIASES
        % Data
        Y_orig
        trueLabels
        y
        % cv_partition
        trIdx
        teIdx
        nfolds
        fold_modeloptions
        hyperparams_by_fold
        % fitted_values
        yfit
        dist_from_hyperplane_xval
        class_probability_xval
        scores_within_id
        Y_within_id
        scorediff
        high_vs_low_scores_within_id
        subjfit
        predictions
        % weights
        w
        mean_vox_weights
        vox_weights
        VOXWEIGHTS
        my_intercepts
        subjbetas
        weight_obj
        boot_w
        boot_w_ste
        boot_w_mean
        wZ
        wP
        wP_fdr_thr
        boot_w_fdrsig
        w_thresh_fdr
        % error_metrics (unwrap .value)
        crossval_accuracy
        crossval_accuracy_within
        d_singleinterval
        d_within
        classification_d_singleinterval
        classification_d_within
        regression_d_singleinterval
        regression_d_within
        prediction_outcome_r
        pred_outcome_r
        cverr
        mse
        rmse
        meanabserr
        r_squared
        var_full
        var_null
        var_reduction
        pred_err
        pred_err_null
        devs_from_full_model
        devs_from_mean_only_model
        r_each_subject
        accuracy
        overallAccuracy
        err
        phi
        % error_metrics descrip (unwrap .descrip)
        crossval_accuracy_descrip
        prediction_outcome_r_descrip
        pred_err_descrip
        var_reduction_descrip
        classification_d_within_descrip
        % stand-alone descrip strings (not paired with a metric)
        accfun_descrip
        cverrfun
        cvoutput_descrip
        error_type
        function_call
        function_handle
        algorithm_name
        % notes
        r_each_subject_note
        r_each_subject_note2
        % ml_model aliases
        ClassificationModel
        SVMModel
        SVRModel
        % bootstrap_results
        WTS
        % diagnostics
        mult_obs_within_person
        ROC_forced_choice
        ROC_single_interval
        all_reg_hyperparams
        % legacy_extras composites (brain wrappers)
        full_model
        data
        covs
    end


    methods

        % -----------------------------------------------------------------
        function obj = predictive_model(varargin)
            % predictive_model  Construct from a struct or from name/value hyperparameters.

            if nargin == 0, return; end

            % --- struct input (legacy wrappers) ---
            if isstruct(varargin{1})
                doverbose = ~any(strcmpi(varargin(2:end), 'noverbose'));
                obj = obj.populate_from_struct(varargin{1}, doverbose);
                return
            end

            % --- name/value hyperparameters (new API) ---
            if ischar(varargin{1}) || isstring(varargin{1})
                p = inputParser;
                p.KeepUnmatched = true;
                addParameter(p, 'algorithm',    obj.algorithm);
                addParameter(p, 'task',         obj.task);
                addParameter(p, 'modeloptions', obj.modeloptions);
                addParameter(p, 'random_state', obj.random_state);
                addParameter(p, 'standardize',  obj.standardize);
                addParameter(p, 'use_parallel', obj.use_parallel);
                addParameter(p, 'cv',           obj.cv);
                addParameter(p, 'scorer',       obj.scorer);
                addParameter(p, 'nboot',        obj.nboot);
                addParameter(p, 'nperm',        obj.nperm);
                addParameter(p, 'do_calibrate', obj.do_calibrate);
                addParameter(p, 'Y_name',       obj.Y_name);
                addParameter(p, 'X_name',       obj.X_name);
                addParameter(p, 'class_labels', obj.class_labels);
                parse(p, varargin{:});
                hp_names = fieldnames(p.Results);
                for i = 1:numel(hp_names)
                    obj.(hp_names{i}) = p.Results.(hp_names{i});
                end
                if ~isempty(fieldnames(p.Unmatched))
                    warning('predictive_model:UnknownHyperparameter', ...
                        'Ignored unknown name/value pairs: %s.', ...
                        strjoin(fieldnames(p.Unmatched), ', '));
                end
                return
            end

            error('predictive_model:InvalidInput', ...
                'First argument must be a struct or a hyperparameter name (char/string).');
        end


        % -----------------------------------------------------------------
        function tf = is_fitted(obj)
            tf = ~isempty(obj.ml_model) ...
                || ~isempty(obj.fold_models) ...
                || ~isempty(fieldnames(obj.fitted_values)) ...
                || ~isempty(fieldnames(obj.weights)) ...
                || ~isempty(fieldnames(obj.error_metrics)) ...
                || ~isempty(fieldnames(obj.cross_classify)) ...
                || ~isempty(fieldnames(obj.bootstrap_results));
        end


        % -----------------------------------------------------------------
        function tf = is_classifier(obj)
            tf = strcmpi(obj.task, 'classification');
            if ~tf && ~isempty(obj.Y)
                tf = numel(unique(obj.Y(~isnan(obj.Y)))) <= 2;
            end
        end


        % -----------------------------------------------------------------
        function tf = is_regressor(obj)
            tf = strcmpi(obj.task, 'regression');
            if ~tf && ~isempty(obj.Y)
                tf = numel(unique(obj.Y(~isnan(obj.Y)))) > 2;
            end
        end


        % -----------------------------------------------------------------
        function new_obj = clone(obj)
            new_obj = predictive_model();
            hp = {'algorithm','task','modeloptions','random_state', ...
                  'standardize','use_parallel','cv','scorer','nboot', ...
                  'nperm','do_calibrate','Y_name','X_name','class_labels'};
            for i = 1:numel(hp)
                new_obj.(hp{i}) = obj.(hp{i});
            end
        end


        % -----------------------------------------------------------------
        function obj = validate_object(obj, varargin)
            doverbose = ~any(strcmpi(varargin, 'noverbose'));
            if ~isempty(obj.algorithm) && ~(ischar(obj.algorithm) || isstring(obj.algorithm))
                error('predictive_model:InvalidProperty', 'algorithm must be char/string.');
            end
            if ~isempty(obj.task) && ~any(strcmpi(obj.task, {'classification','regression'}))
                error('predictive_model:InvalidProperty', ...
                    'task must be ''classification'' or ''regression''.');
            end
            if ~isempty(obj.fit_type) && ~any(strcmpi(obj.fit_type, {'crossval','insample','test'}))
                error('predictive_model:InvalidProperty', ...
                    'fit_type must be ''crossval'', ''insample'', or ''test''.');
            end
            if ~isempty(obj.modeloptions)
                if ~iscell(obj.modeloptions)
                    error('predictive_model:InvalidProperty', ...
                        'modeloptions must be a cell array of name/value pairs.');
                end
                names = obj.modeloptions(1:2:end);
                bad = ~cellfun(@(c) ischar(c) || isstring(c), names);
                if any(bad)
                    error('predictive_model:InvalidProperty', ...
                        'modeloptions: option names (odd positions) must be char/string.');
                end
            end
            cats = {'inputParameters','cv_partition','fitted_values','weights', ...
                    'error_metrics','bootstrap_results','permutation_results', ...
                    'diagnostics','cross_classify','legacy_extras'};
            for i = 1:numel(cats)
                if ~isstruct(obj.(cats{i}))
                    error('predictive_model:InvalidProperty', '%s must be a struct.', cats{i});
                end
            end
            if ~isempty(obj.Y)
                % Y may be a numeric vector (single-dataset wrappers) OR
                % a cell array of numeric vectors (multi-subject wrappers).
                if ~(isnumeric(obj.Y) || iscell(obj.Y))
                    error('predictive_model:InvalidProperty', ...
                        'Y must be a numeric vector or a cell array of numeric vectors.');
                end
            end
            if doverbose, disp('predictive_model object validated successfully.'); end
        end


        % -----------------------------------------------------------------
        % LEGACY DEPENDENT ALIASES (read-only)
        % -----------------------------------------------------------------
        % Data
        function v = get.Y_orig(obj),     v = obj.Y; end
        function v = get.y(obj),          v = obj.Y; end
        function v = get.trueLabels(obj), v = predictive_model.field_or_empty(obj.fitted_values, 'Y_per_fold'); end

        % cv_partition
        function v = get.trIdx(obj),               v = predictive_model.field_or_empty(obj.cv_partition, 'trIdx'); end
        function v = get.teIdx(obj),               v = predictive_model.field_or_empty(obj.cv_partition, 'teIdx'); end
        function v = get.nfolds(obj),              v = predictive_model.field_or_empty(obj.cv_partition, 'nfolds'); end
        function v = get.fold_modeloptions(obj),   v = predictive_model.field_or_empty(obj.cv_partition, 'fold_modeloptions'); end
        function v = get.hyperparams_by_fold(obj), v = predictive_model.field_or_empty(obj.cv_partition, 'hyperparams_by_fold'); end

        % fitted_values
        function v = get.yfit(obj),                         v = predictive_model.field_or_empty(obj.fitted_values, 'yfit'); end
        function v = get.dist_from_hyperplane_xval(obj)
            % Returns scores when score_type='distance', else legacy field.
            if isfield(obj.fitted_values, 'score_type') && strcmpi(obj.fitted_values.score_type, 'distance')
                v = predictive_model.field_or_empty(obj.fitted_values, 'scores');
            else
                v = predictive_model.field_or_empty(obj.fitted_values, 'dist_from_hyperplane_xval');
            end
        end
        function v = get.class_probability_xval(obj)
            if isfield(obj.fitted_values, 'score_type') && strcmpi(obj.fitted_values.score_type, 'probability')
                v = predictive_model.field_or_empty(obj.fitted_values, 'scores');
            else
                v = predictive_model.field_or_empty(obj.fitted_values, 'class_probability_xval');
            end
        end
        function v = get.scores_within_id(obj),             v = predictive_model.field_or_empty(obj.fitted_values, 'scores_within_id'); end
        function v = get.Y_within_id(obj),                  v = predictive_model.field_or_empty(obj.fitted_values, 'Y_within_id'); end
        function v = get.scorediff(obj),                    v = predictive_model.field_or_empty(obj.fitted_values, 'scorediff'); end
        function v = get.high_vs_low_scores_within_id(obj), v = predictive_model.field_or_empty(obj.fitted_values, 'high_vs_low_scores_within_id'); end
        function v = get.subjfit(obj),                      v = predictive_model.field_or_empty(obj.fitted_values, 'subjfit'); end
        function v = get.predictions(obj),                  v = predictive_model.field_or_empty(obj.fitted_values, 'predictions'); end

        % weights (and merged former weight_stats)
        function v = get.w(obj),                v = predictive_model.field_or_empty(obj.weights, 'w'); end
        function v = get.mean_vox_weights(obj), v = predictive_model.field_or_empty(obj.weights, 'mean_vox_weights'); end
        function v = get.vox_weights(obj),      v = predictive_model.field_or_empty(obj.weights, 'w_perfold'); end
        function v = get.VOXWEIGHTS(obj),       v = predictive_model.field_or_empty(obj.weights, 'w_bootstrap'); end
        function v = get.my_intercepts(obj),    v = predictive_model.field_or_empty(obj.weights, 'my_intercepts'); end
        function v = get.subjbetas(obj),        v = predictive_model.field_or_empty(obj.weights, 'subjbetas'); end
        function v = get.weight_obj(obj),       v = predictive_model.field_or_empty(obj.weights, 'weight_obj'); end
        function v = get.boot_w(obj),           v = predictive_model.field_or_empty(obj.weights, 'boot_w'); end
        function v = get.boot_w_ste(obj),       v = predictive_model.field_or_empty(obj.weights, 'boot_w_ste'); end
        function v = get.boot_w_mean(obj),      v = predictive_model.field_or_empty(obj.weights, 'boot_w_mean'); end
        function v = get.wZ(obj),               v = predictive_model.field_or_empty(obj.weights, 'z'); end
        function v = get.wP(obj),               v = predictive_model.field_or_empty(obj.weights, 'p'); end
        function v = get.wP_fdr_thr(obj),       v = predictive_model.field_or_empty(obj.weights, 'fdr_thr'); end
        function v = get.boot_w_fdrsig(obj),    v = predictive_model.field_or_empty(obj.weights, 'fdr_sig'); end
        function v = get.w_thresh_fdr(obj),     v = predictive_model.field_or_empty(obj.weights, 'thresh_fdr'); end

        % error_metrics — unwrap .value for arithmetic
        function v = get.crossval_accuracy(obj),                v = predictive_model.metric_value(obj.error_metrics, 'crossval_accuracy'); end
        function v = get.crossval_accuracy_within(obj),         v = predictive_model.metric_value(obj.error_metrics, 'crossval_accuracy_within'); end
        function v = get.d_singleinterval(obj),                 v = predictive_model.metric_value(obj.error_metrics, 'd_singleinterval'); end
        function v = get.d_within(obj),                         v = predictive_model.metric_value(obj.error_metrics, 'd_within'); end
        function v = get.classification_d_singleinterval(obj),  v = predictive_model.metric_value(obj.error_metrics, 'd_singleinterval'); end
        function v = get.classification_d_within(obj),          v = predictive_model.metric_value(obj.error_metrics, 'd_within'); end
        function v = get.regression_d_singleinterval(obj),      v = predictive_model.metric_value(obj.error_metrics, 'd_singleinterval'); end
        function v = get.regression_d_within(obj),              v = predictive_model.metric_value(obj.error_metrics, 'd_within'); end
        function v = get.prediction_outcome_r(obj),             v = predictive_model.metric_value(obj.error_metrics, 'prediction_outcome_r'); end
        function v = get.pred_outcome_r(obj),                   v = predictive_model.metric_value(obj.error_metrics, 'prediction_outcome_r'); end
        function v = get.cverr(obj),                            v = predictive_model.metric_value(obj.error_metrics, 'cverr'); end
        function v = get.mse(obj),                              v = predictive_model.metric_value(obj.error_metrics, 'mse'); end
        function v = get.rmse(obj),                             v = predictive_model.metric_value(obj.error_metrics, 'rmse'); end
        function v = get.meanabserr(obj),                       v = predictive_model.metric_value(obj.error_metrics, 'meanabserr'); end
        function v = get.r_squared(obj),                        v = predictive_model.metric_value(obj.error_metrics, 'r_squared'); end
        function v = get.var_full(obj),                         v = predictive_model.metric_value(obj.error_metrics, 'var_full'); end
        function v = get.var_null(obj),                         v = predictive_model.metric_value(obj.error_metrics, 'var_null'); end
        function v = get.var_reduction(obj),                    v = predictive_model.metric_value(obj.error_metrics, 'var_reduction'); end
        function v = get.pred_err(obj),                         v = predictive_model.metric_value(obj.error_metrics, 'pred_err'); end
        function v = get.pred_err_null(obj),                    v = predictive_model.metric_value(obj.error_metrics, 'pred_err_null'); end
        function v = get.devs_from_full_model(obj),             v = predictive_model.metric_value(obj.error_metrics, 'devs_from_full_model'); end
        function v = get.devs_from_mean_only_model(obj),        v = predictive_model.metric_value(obj.error_metrics, 'devs_from_mean_only_model'); end
        function v = get.r_each_subject(obj),                   v = predictive_model.metric_value(obj.error_metrics, 'r_each_subject'); end
        function v = get.accuracy(obj),                         v = predictive_model.metric_value(obj.error_metrics, 'accuracy'); end
        function v = get.overallAccuracy(obj),                  v = predictive_model.metric_value(obj.error_metrics, 'overallAccuracy'); end
        function v = get.err(obj),                              v = predictive_model.metric_value(obj.error_metrics, 'err'); end
        function v = get.phi(obj),                              v = predictive_model.metric_value(obj.error_metrics, 'phi'); end

        % error_metrics — descrip side of the tuple
        function v = get.crossval_accuracy_descrip(obj),        v = predictive_model.metric_descrip(obj.error_metrics, 'crossval_accuracy'); end
        function v = get.prediction_outcome_r_descrip(obj),     v = predictive_model.metric_descrip(obj.error_metrics, 'prediction_outcome_r'); end
        function v = get.pred_err_descrip(obj),                 v = predictive_model.metric_descrip(obj.error_metrics, 'pred_err'); end
        function v = get.var_reduction_descrip(obj),            v = predictive_model.metric_descrip(obj.error_metrics, 'var_reduction'); end
        function v = get.classification_d_within_descrip(obj),  v = predictive_model.metric_descrip(obj.error_metrics, 'd_within'); end

        % notes
        function v = get.r_each_subject_note(obj)
            v = predictive_model.lookup_note(obj.note, 'r_each_subject');
        end
        function v = get.r_each_subject_note2(obj)
            v = predictive_model.lookup_note(obj.note, 'r_each_subject_2');
        end

        % stand-alone descriptions (not paired with a metric)
        function v = get.accfun_descrip(obj),                   v = predictive_model.field_or_empty(obj.legacy_extras, 'accfun_descrip'); end
        function v = get.cverrfun(obj),                         v = predictive_model.field_or_empty(obj.legacy_extras, 'cverrfun'); end
        function v = get.cvoutput_descrip(obj),                 v = predictive_model.field_or_empty(obj.legacy_extras, 'cvoutput_descrip'); end
        function v = get.error_type(obj),                       v = predictive_model.field_or_empty(obj.legacy_extras, 'error_type'); end
        function v = get.function_call(obj),                    v = predictive_model.field_or_empty(obj.legacy_extras, 'function_call'); end
        function v = get.function_handle(obj),                  v = predictive_model.field_or_empty(obj.legacy_extras, 'function_handle'); end
        function v = get.algorithm_name(obj),                   v = predictive_model.field_or_empty(obj.legacy_extras, 'algorithm_name'); end

        % ml_model aliases
        function v = get.ClassificationModel(obj), v = obj.ml_model; end
        function v = get.SVMModel(obj),            v = obj.ml_model; end
        function v = get.SVRModel(obj),            v = obj.ml_model; end

        % bootstrap_results
        function v = get.WTS(obj), v = predictive_model.field_or_empty(obj.bootstrap_results, 'WTS'); end

        % diagnostics
        function v = get.mult_obs_within_person(obj), v = predictive_model.field_or_empty(obj.diagnostics, 'mult_obs_within_person'); end
        function v = get.ROC_forced_choice(obj),      v = predictive_model.field_or_empty(obj.diagnostics, 'ROC_forced_choice'); end
        function v = get.ROC_single_interval(obj),    v = predictive_model.field_or_empty(obj.diagnostics, 'ROC_single_interval'); end
        function v = get.all_reg_hyperparams(obj),    v = predictive_model.field_or_empty(obj.diagnostics, 'all_reg_hyperparams'); end

        % legacy_extras composite (brain wrappers)
        function v = get.full_model(obj), v = predictive_model.field_or_empty(obj.legacy_extras, 'full_model'); end
        function v = get.data(obj),       v = predictive_model.field_or_empty(obj.legacy_extras, 'data'); end
        function v = get.covs(obj),       v = predictive_model.field_or_empty(obj.legacy_extras, 'covs'); end

    end


    methods (Access = protected)

        % -----------------------------------------------------------------
        function obj = populate_from_struct(obj, S, doverbose)
            % Walk input struct, route each field to its categorised home.

            routing = predictive_model.field_routing();
            hp_props = predictive_model.hyperparameter_names();
            keys = routing(:, 1);

            fn = fieldnames(S);
            unrouted_to_extras = {};

            for i = 1:numel(fn)
                name = fn{i};
                val  = S.(name);

                % 1) hyperparameter direct match
                if ismember(name, hp_props)
                    obj.(name) = val;
                    continue
                end

                % 2) top-level new properties + a few legacy top-level
                if any(strcmp(name, {'Y','id','omitted_cases','omitted_features', ...
                                     'fit_type','accfun','history','removed_voxels'}))
                    obj.(name) = val;
                    continue
                end

                % 3) routing table (dotted path target)
                idx = find(strcmp(keys, name), 1);
                if ~isempty(idx)
                    target = routing{idx, 2};
                    obj = predictive_model.write_path(obj, target, val);
                    continue
                end

                % 4) catch-all
                obj.legacy_extras.(name) = val;
                unrouted_to_extras{end+1, 1} = name; %#ok<AGROW>
            end

            obj = obj.validate_object('noverbose');

            if doverbose && ~isempty(unrouted_to_extras)
                fprintf(['predictive_model: %d input field(s) had no categorised home ' ...
                         'and were stored in obj.legacy_extras:\n'], ...
                        numel(unrouted_to_extras));
                disp(unrouted_to_extras);
            end
        end

    end


    methods (Static)

        % -----------------------------------------------------------------
        function v = field_or_empty(s, name)
            if isstruct(s) && isfield(s, name)
                v = s.(name);
            else
                v = [];
            end
        end


        % -----------------------------------------------------------------
        function v = metric_value(em_struct, name)
            % Unwrap .value from a (value, descrip) tuple.
            if isstruct(em_struct) && isfield(em_struct, name)
                entry = em_struct.(name);
                if isstruct(entry) && isfield(entry, 'value')
                    v = entry.value;
                else
                    v = entry;
                end
            else
                v = [];
            end
        end


        % -----------------------------------------------------------------
        function v = metric_descrip(em_struct, name)
            % Unwrap .descrip from a (value, descrip) tuple.
            if isstruct(em_struct) && isfield(em_struct, name)
                entry = em_struct.(name);
                if isstruct(entry) && isfield(entry, 'descrip')
                    v = entry.descrip;
                else
                    v = '';
                end
            else
                v = '';
            end
        end


        % -----------------------------------------------------------------
        function names = hyperparameter_names()
            names = {'algorithm','task','modeloptions','random_state', ...
                     'standardize','use_parallel','cv','scorer','nboot', ...
                     'nperm','do_calibrate','Y_name','X_name','class_labels'};
        end


        % -----------------------------------------------------------------
        function v = lookup_note(note_cell, prefix)
            % Find a note string starting with `prefix:` and return its body.
            v = '';
            if ~iscell(note_cell), return; end
            for i = 1:numel(note_cell)
                s = note_cell{i};
                if ischar(s) || isstring(s)
                    s = char(s);
                    if startsWith(s, [prefix ':'])
                        v = strtrim(s(numel(prefix)+2:end));
                        return
                    end
                end
            end
        end


        % -----------------------------------------------------------------
        function obj = write_path(obj, target_path, value)
            % Write `value` into obj at a dotted path like
            % 'error_metrics.pred_err.value' or 'weights.z'.
            %
            % An embedded '+' marks a cell-array append:
            %   'note+'           append the value as a new cell entry
            %   'note+r_each_subject'  append with prefix "r_each_subject: "

            if contains(target_path, '+')
                tokens = strsplit(target_path, '+');
                target_path = tokens{1};
                prefix = tokens{2};
                if ~ischar(value) && ~isstring(value)
                    value = char(string(value));
                end
                value = char(value);
                if ~isempty(prefix)
                    value = [prefix ': ' value];
                end
                parts = strsplit(target_path, '.');
                top = parts{1};
                if numel(parts) == 1
                    cur = obj.(top);
                    if ~iscell(cur), cur = {}; end
                    cur{end+1, 1} = value;
                    obj.(top) = cur;
                    return
                end
            end

            parts = strsplit(target_path, '.');
            top = parts{1};

            if numel(parts) == 1
                obj.(top) = value;
                return
            end

            if ~isstruct(obj.(top))
                obj.(top) = struct();
            end
            sub_struct = obj.(top);
            sub_struct = predictive_model.setfield_recursive(sub_struct, parts(2:end), value);
            obj.(top) = sub_struct;
        end


        % -----------------------------------------------------------------
        function s = setfield_recursive(s, parts, value)
            if numel(parts) == 1
                s.(parts{1}) = value;
            else
                if ~isfield(s, parts{1}) || ~isstruct(s.(parts{1}))
                    s.(parts{1}) = struct();
                end
                s.(parts{1}) = predictive_model.setfield_recursive(s.(parts{1}), parts(2:end), value);
            end
        end


        % -----------------------------------------------------------------
        function [omitted_cases, omitted_features] = detect_bad_data(X, Y)
            % detect_bad_data  Pre-fit data-quality check.
            %
            % Returns logical vectors over the ORIGINAL cases and features:
            %   omitted_cases     true for rows with any NaN/Inf in X or Y
            %   omitted_features  true for columns that, after dropping bad
            %                     cases, are all-NaN, constant (zero variance),
            %                     or contain Inf.
            %
            % Wrappers apply these by:
            %     X(omitted_cases, :) = []; Y(omitted_cases) = [];
            %     X(:, omitted_features) = [];
            %     pm.omitted_cases    = omitted_cases;
            %     pm.omitted_features = omitted_features;

            Y = Y(:);
            omitted_cases = isnan(Y) | isinf(Y) | any(isnan(X), 2) | any(isinf(X), 2);

            Xg = X(~omitted_cases, :);
            if isempty(Xg)
                omitted_features = false(size(X, 2), 1);
                return
            end
            all_nan = all(isnan(Xg), 1);
            any_inf = any(isinf(Xg), 1);
            zero_var = false(1, size(Xg, 2));
            ok_cols = ~all_nan;
            if any(ok_cols)
                zero_var(ok_cols) = var(Xg(:, ok_cols), 0, 1, 'omitnan') == 0;
            end
            omitted_features = (all_nan | any_inf | zero_var).';
        end


        % -----------------------------------------------------------------
        function routing = field_routing()
            % field_routing  Lookup table: input-struct field name -> dotted
            % target path inside the predictive_model object.
            %
            % An optional trailing '+' on a path means "append to cell array"
            % (used for `note`).

            routing = { ...
                % --- DATA (top-level after consolidation) ---
                'Y',                                 'Y'
                'y',                                 'Y'
                'Y_orig',                            'Y'
                'id',                                'id'
                'omitted_cases',                     'omitted_cases'
                'omitted_features',                  'omitted_features'

                % --- inputParameters (only canonical name accepted; INPUTS / inputs / inputOptions deprecated) ---
                'inputParameters',                   'inputParameters'

                % --- fit metadata ---
                'fit_type',                          'fit_type'

                % --- cv_partition ---
                'trIdx',                             'cv_partition.trIdx'
                'teIdx',                             'cv_partition.teIdx'
                'nfolds',                            'cv_partition.nfolds'
                'fold_modeloptions',                 'cv_partition.fold_modeloptions'
                'hyperparams_by_fold',               'cv_partition.hyperparams_by_fold'
                'cvpartition',                       'cv_partition.cvpartition'

                % --- fitted_values ---
                'yfit',                              'fitted_values.yfit'
                'scores',                            'fitted_values.scores'
                'score_type',                        'fitted_values.score_type'
                'dist_from_hyperplane_xval',         'fitted_values.dist_from_hyperplane_xval'
                'class_probability_xval',            'fitted_values.class_probability_xval'
                'scores_within_id',                  'fitted_values.scores_within_id'
                'Y_within_id',                       'fitted_values.Y_within_id'
                'scorediff',                         'fitted_values.scorediff'
                'high_vs_low_scores_within_id',      'fitted_values.high_vs_low_scores_within_id'
                'subjfit',                           'fitted_values.subjfit'
                'predictions',                       'fitted_values.predictions'
                'trueLabels',                        'fitted_values.Y_per_fold'

                % --- weights (incl. former weight_stats and renamed VOXWEIGHTS/vox_weights) ---
                'w',                                 'weights.w'
                'mean_vox_weights',                  'weights.mean_vox_weights'
                'vox_weights',                       'weights.w_perfold'
                'VOXWEIGHTS',                        'weights.w_bootstrap'
                'my_intercepts',                     'weights.my_intercepts'
                'subjbetas',                         'weights.subjbetas'
                'weight_obj',                        'weights.weight_obj'
                'boot_w',                            'weights.boot_w'
                'boot_w_ste',                        'weights.boot_w_ste'
                'boot_w_mean',                       'weights.boot_w_mean'
                'wZ',                                'weights.z'
                'wP',                                'weights.p'
                'wP_fdr_thr',                        'weights.fdr_thr'
                'boot_w_fdrsig',                     'weights.fdr_sig'
                'w_thresh_fdr',                      'weights.thresh_fdr'

                % --- error_metrics — value half of (value, descrip) tuple ---
                'crossval_accuracy',                 'error_metrics.crossval_accuracy.value'
                'crossval_accuracy_within',          'error_metrics.crossval_accuracy_within.value'
                'd_singleinterval',                  'error_metrics.d_singleinterval.value'
                'd_within',                          'error_metrics.d_within.value'
                'classification_d_singleinterval',   'error_metrics.d_singleinterval.value'
                'regression_d_singleinterval',       'error_metrics.d_singleinterval.value'
                'classification_d_within',           'error_metrics.d_within.value'
                'regression_d_within',               'error_metrics.d_within.value'
                'prediction_outcome_r',              'error_metrics.prediction_outcome_r.value'
                'pred_outcome_r',                    'error_metrics.prediction_outcome_r.value'
                'cverr',                             'error_metrics.cverr.value'
                'mse',                               'error_metrics.mse.value'
                'rmse',                              'error_metrics.rmse.value'
                'meanabserr',                        'error_metrics.meanabserr.value'
                'r_squared',                         'error_metrics.r_squared.value'
                'var_full',                          'error_metrics.var_full.value'
                'var_null',                          'error_metrics.var_null.value'
                'var_reduction',                     'error_metrics.var_reduction.value'
                'pred_err',                          'error_metrics.pred_err.value'
                'pred_err_null',                     'error_metrics.pred_err_null.value'
                'devs_from_full_model',              'error_metrics.devs_from_full_model.value'
                'devs_from_mean_only_model',         'error_metrics.devs_from_mean_only_model.value'
                'r_each_subject',                    'error_metrics.r_each_subject.value'
                'accuracy',                          'error_metrics.accuracy.value'
                'overallAccuracy',                   'error_metrics.overallAccuracy.value'
                'err',                               'error_metrics.err.value'
                'phi',                               'error_metrics.phi.value'

                % --- error_metrics — descrip half ---
                'crossval_accuracy_descrip',         'error_metrics.crossval_accuracy.descrip'
                'prediction_outcome_r_descrip',      'error_metrics.prediction_outcome_r.descrip'
                'pred_err_descrip',                  'error_metrics.pred_err.descrip'
                'var_reduction_descrip',             'error_metrics.var_reduction.descrip'
                'classification_d_within_descrip',   'error_metrics.d_within.descrip'

                % --- notes (cell-array append) ---
                'note',                              'note+'
                'r_each_subject_note',               'note+r_each_subject'
                'r_each_subject_note2',              'note+r_each_subject_2'

                % --- ml_model (legacy single-object aliases) ---
                'ClassificationModel',               'ml_model'
                'SVMModel',                          'ml_model'
                'SVRModel',                          'ml_model'

                % --- fold_models ---
                'models',                            'fold_models'

                % --- bootstrap_results ---
                'bootstrap',                         'bootstrap_results'
                'WTS',                               'bootstrap_results.WTS'
                'boot_weights',                      'bootstrap_results.boot_weights'

                % --- diagnostics ---
                'ROC_forced_choice',                 'diagnostics.ROC_forced_choice'
                'ROC_single_interval',               'diagnostics.ROC_single_interval'
                'mult_obs_within_person',            'diagnostics.mult_obs_within_person'
                'all_reg_hyperparams',               'diagnostics.all_reg_hyperparams'

                % --- cross_classify ---
                'stats1',                            'cross_classify.stats1'
                'stats2',                            'cross_classify.stats2'
                'test_results',                      'cross_classify.test_results'
                'crosstestfit',                      'cross_classify.crosstestfit'
                'cvoutput',                          'cross_classify.cvoutput'
                'test_Y',                            'cross_classify.test_Y'
                'all',                               'cross_classify.all'

                % --- stand-alone descrip strings (kept in legacy_extras) ---
                'accfun_descrip',                    'legacy_extras.accfun_descrip'
                'cverrfun',                          'legacy_extras.cverrfun'
                'cvoutput_descrip',                  'legacy_extras.cvoutput_descrip'
                'error_type',                        'legacy_extras.error_type'
                'function_call',                     'legacy_extras.function_call'
                'function_handle',                   'legacy_extras.function_handle'
                'algorithm_name',                    'legacy_extras.algorithm_name'

                % --- legacy_extras catch-all for brain-wrapper composite fields ---
                'covs',                              'legacy_extras.covs'
                'data',                              'legacy_extras.data'
                'full_model',                        'legacy_extras.full_model'
            };
        end

    end

end
