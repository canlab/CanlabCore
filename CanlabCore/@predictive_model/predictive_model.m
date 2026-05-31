classdef predictive_model
    % predictive_model  Object representing a multivariate predictive model.
    %
    % This class is the canonical output of CanlabCore's predictive-modeling
    % entry points (xval_SVM, xval_SVR, xval_discriminant_classifier,
    % xval_cross_classify, xval_regression_multisubject, xval_lasso_brain,
    % xval_ridge_brain, xval_bestsubsets_brain, fmri_data.predict, ...).
    %
    % Properties are organised into two top-level blocks:
    %
    %   HYPERPARAMETERS   (set by the user / wrapper)
    %     algorithm, task, modeloptions, random_state, standardize,
    %     use_parallel, cv, scorer, nboot, nperm, do_calibrate,
    %     Y_name, X_name, class_labels.
    %
    %   FITTED STATE      (populated by fit / crossval / bootstrap / ...)
    %     The fitted-state surface is grouped into 13 named sub-struct
    %     properties plus a few top-level fields, so the class stays
    %     grok-able as more algorithm families are added.
    %
    %       inputs               .Y .id .Y_orig .INPUTS
    %       cv_partition         .trIdx .teIdx .nfolds .fold_modeloptions
    %                            .hyperparams_by_fold
    %       fitted_values        .yfit .dist_from_hyperplane_xval
    %                            .class_probability_xval .scores_within_id
    %                            .Y_within_id .scorediff .subjfit
    %                            .high_vs_low_scores_within_id .predictions
    %                            .trueLabels
    %       weights              .w (= Beta) .mean_vox_weights .vox_weights
    %                            .VOXWEIGHTS .my_intercepts .subjbetas
    %                            .weight_obj   (statistic_image when available)
    %       weight_stats         .boot_w .boot_w_ste .boot_w_mean .wZ .wP
    %                            .wP_fdr_thr .boot_w_fdrsig .w_thresh_fdr
    %       error_metrics        .crossval_accuracy .crossval_accuracy_within
    %                            .d_singleinterval .d_within
    %                            .prediction_outcome_r .cverr .mse .rmse
    %                            .meanabserr .r_squared .var_full .var_null
    %                            .var_reduction .pred_err .pred_err_null
    %                            .devs_from_full_model
    %                            .devs_from_mean_only_model .r_each_subject
    %                            .accuracy .overallAccuracy .err
    %       descriptions         .accfun_descrip .cverrfun
    %                            .crossval_accuracy_descrip
    %                            .prediction_outcome_r_descrip
    %                            .pred_err_descrip .var_reduction_descrip
    %                            .classification_d_within_descrip
    %                            .note .r_each_subject_note
    %                            .r_each_subject_note2 .cvoutput_descrip
    %                            .error_type .function_call .function_handle
    %                            .algorithm_name
    %       ml_model             trained MATLAB model object (e.g.,
    %                            ClassificationSVM, RegressionLinear, ...);
    %                            accepts legacy field names SVMModel,
    %                            SVRModel, ClassificationModel.
    %       fold_models          cell array of per-fold trained MATLAB models
    %                            (accepts legacy field name "models").
    %       bootstrap_results    full bootstrap output (e.g., the .bootstrap
    %                            sub-struct produced by xval_lasso_brain);
    %                            .WTS is also routed here.
    %       permutation_results  full permutation-test output
    %       diagnostics          .ROC_forced_choice .ROC_single_interval
    %                            .mult_obs_within_person
    %                            .all_reg_hyperparams .phi
    %       cross_classify       cross-classification output from
    %                            xval_cross_classify (.stats1 .stats2
    %                            .test_results .crosstestfit .cvoutput
    %                            .test_Y .all)
    %       legacy_extras        catch-all struct for any input-struct field
    %                            that does not yet have a categorised home;
    %                            stored verbatim so nothing is lost.
    %       accfun               objective function handle used in CV
    %       history              cell of one-line strings, image_vector-style
    %       removed_voxels       logical vector, image_vector-style
    %
    %   LEGACY ALIASES (read-only, Dependent)
    %     For backward compatibility with code written against the flat
    %     output struct that xval_SVM / xval_SVR / etc. used to return, the
    %     class exposes Dependent properties at the top level that forward
    %     to the categorised sub-struct fields:
    %
    %       Y, id, yfit, w, boot_w, boot_w_ste, boot_w_mean, wZ, wP,
    %       wP_fdr_thr, boot_w_fdrsig, w_thresh_fdr,
    %       dist_from_hyperplane_xval, class_probability_xval,
    %       crossval_accuracy, crossval_accuracy_within,
    %       classification_d_singleinterval, classification_d_within,
    %       regression_d_singleinterval, regression_d_within,
    %       d_singleinterval, d_within, mult_obs_within_person,
    %       ClassificationModel, SVMModel, SVRModel, trIdx, teIdx, nfolds,
    %       accfun_descrip, cverrfun, cverr, crossval_accuracy_descrip,
    %       prediction_outcome_r, prediction_outcome_r_descrip,
    %       pred_outcome_r, mse, rmse, meanabserr, r_squared,
    %       scorediff, scores_within_id, Y_within_id,
    %       ROC_forced_choice, ROC_single_interval, weight_obj, modeloptions.
    %
    %     Reading these (`pm.w`, `pm.yfit`, ...) keeps working exactly as
    %     before. Writing to them is not supported; mutate the categorised
    %     sub-structs directly via methods.
    %
    % :Usage:
    % ::
    %     pm = predictive_model();             % empty object
    %     pm = predictive_model(S);            % from struct (legacy wrappers)
    %     pm = predictive_model(S, 'noverbose'); % suppress "fields not copied"
    %     pm = predictive_model('algorithm','svm','task','classification', ...
    %                           'modeloptions',{'KernelFunction','linear'});
    %
    % When constructed from a struct (the path used by every xval_*.m wrapper
    % today), the constructor walks the struct's fields and routes each one
    % into its categorised home. Legacy translations:
    %     SVMModel, SVRModel                 -> ml_model
    %     ClassificationModel                -> ml_model
    %     classification_d_singleinterval    -> error_metrics.d_singleinterval
    %     regression_d_singleinterval        -> error_metrics.d_singleinterval
    %     classification_d_within            -> error_metrics.d_within
    %     regression_d_within                -> error_metrics.d_within
    %     pred_outcome_r                     -> error_metrics.prediction_outcome_r
    %     models                             -> fold_models
    %     bootstrap                          -> bootstrap_results
    %
    % :Examples:
    % ::
    %     % Build from a struct produced by xval_SVM:
    %     S  = xval_SVM(X, Y, id, 'nooptimize', 'norepeats', 'nobootstrap');
    %     pm = predictive_model(S);
    %     pm.is_fitted               % true
    %     pm.is_classifier           % true (Y has two unique values)
    %     pm.crossval_accuracy       % legacy alias -> error_metrics.crossval_accuracy
    %     pm.fitted_values.yfit      % canonical path to held-out predictions
    %     pm.weights.w               % canonical path to model weights
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
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    % -------------------------------------------------------------------------

    properties  % HYPERPARAMETERS
        algorithm       = '';     % string, e.g. 'svm','svr','lassopcr','ridge','lda'
        task            = '';     % 'classification' | 'regression'
        modeloptions    = {};     % cell of name/value pairs forwarded to fit fn
        random_state    = [];     % rng seed for reproducibility
        standardize     = [];     % logical
        use_parallel    = [];     % logical
        cv              = [];     % splitter object (Phase 1)
        scorer          = [];     % scorer object (Phase 1)
        nboot           = [];     % bootstrap samples (default chosen at fit)
        nperm           = [];     % permutation samples
        do_calibrate    = [];     % logical; Platt-style calibration
        Y_name          = '';
        X_name          = '';
        class_labels    = {};
    end

    properties (SetAccess = protected)  % CATEGORISED FITTED STATE
        inputs               = struct();
        cv_partition         = struct();
        fitted_values        = struct();
        weights              = struct();
        weight_stats         = struct();
        error_metrics        = struct();
        descriptions         = struct();
        ml_model             = [];
        fold_models          = {};
        bootstrap_results    = struct();
        permutation_results  = struct();
        diagnostics          = struct();
        cross_classify       = struct();
        legacy_extras        = struct();
        accfun               = [];    % top-level: objective function handle
        history              = {};
        removed_voxels       = [];
    end

    properties (Dependent, SetAccess = private)  % LEGACY READ-ONLY ALIASES
        % -- inputs --
        Y
        id
        Y_orig
        INPUTS
        % -- cv_partition --
        trIdx
        teIdx
        nfolds
        fold_modeloptions
        hyperparams_by_fold
        % -- fitted_values --
        yfit
        dist_from_hyperplane_xval
        class_probability_xval
        scores_within_id
        Y_within_id
        scorediff
        high_vs_low_scores_within_id
        subjfit
        predictions
        trueLabels
        % -- weights --
        w
        mean_vox_weights
        vox_weights
        VOXWEIGHTS
        my_intercepts
        subjbetas
        weight_obj
        % -- weight_stats --
        boot_w
        boot_w_ste
        boot_w_mean
        wZ
        wP
        wP_fdr_thr
        boot_w_fdrsig
        w_thresh_fdr
        % -- error_metrics (incl. legacy d-prefix variants) --
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
        % -- descriptions --
        accfun_descrip
        cverrfun
        crossval_accuracy_descrip
        prediction_outcome_r_descrip
        pred_err_descrip
        var_reduction_descrip
        classification_d_within_descrip
        note
        r_each_subject_note
        r_each_subject_note2
        cvoutput_descrip
        error_type
        function_call
        function_handle
        algorithm_name
        % -- ml_model legacy single-object aliases --
        ClassificationModel
        SVMModel
        SVRModel
        % -- bootstrap_results --
        WTS
        % -- diagnostics --
        mult_obs_within_person
        ROC_forced_choice
        ROC_single_interval
        all_reg_hyperparams
        % -- legacy_extras (brain-wrapper composite fields) --
        full_model
        data
        covs
        inputOptions
    end


    methods

        % -----------------------------------------------------------------
        function obj = predictive_model(varargin)
            % predictive_model  Construct a predictive_model object.
            %
            % Usage:
            %   pm = predictive_model();
            %   pm = predictive_model(S);                   % struct path
            %   pm = predictive_model(S, 'noverbose');
            %   pm = predictive_model('algorithm','svm', 'task','classification', ...);

            if nargin == 0, return; end

            % --- Path 1: struct input (legacy wrappers) ---
            if isstruct(varargin{1})
                doverbose = ~any(strcmpi(varargin(2:end), 'noverbose'));
                obj = obj.populate_from_struct(varargin{1}, doverbose);
                return
            end

            % --- Path 2: name-value pairs (new API) ---
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
                        ['Ignored unknown name/value pairs: %s. ' ...
                         'For struct input, pass the struct as the first argument.'], ...
                        strjoin(fieldnames(p.Unmatched), ', '));
                end
                return
            end

            error('predictive_model:InvalidInput', ...
                'First argument must be a struct or a hyperparameter name (char/string).');
        end


        % -----------------------------------------------------------------
        function tf = is_fitted(obj)
            % is_fitted  True if any categorised fitted-state field has content.
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
            % is_classifier  True if task is 'classification' or Y is binary.
            tf = strcmpi(obj.task, 'classification');
            if ~tf && isfield(obj.inputs, 'Y') && ~isempty(obj.inputs.Y)
                Yv = obj.inputs.Y;
                tf = numel(unique(Yv(~isnan(Yv)))) <= 2;
            end
        end


        % -----------------------------------------------------------------
        function tf = is_regressor(obj)
            % is_regressor  True if task is 'regression' (or Y has >2 unique values).
            tf = strcmpi(obj.task, 'regression');
            if ~tf && isfield(obj.inputs, 'Y') && ~isempty(obj.inputs.Y)
                Yv = obj.inputs.Y;
                tf = numel(unique(Yv(~isnan(Yv)))) > 2;
            end
        end


        % -----------------------------------------------------------------
        function new_obj = clone(obj)
            % clone  Return a fresh predictive_model with the same
            %        hyperparameters but no fitted state.
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
            % validate_object  Light type/shape checks for populated fields.
            %
            % Categorised sub-structs are validated only when non-empty.
            % Pass 'noverbose' to suppress the success message.

            doverbose = ~any(strcmpi(varargin, 'noverbose'));

            % Hyperparameter shape checks (only when set).
            if ~isempty(obj.algorithm) && ~(ischar(obj.algorithm) || isstring(obj.algorithm))
                error('predictive_model:InvalidProperty', 'algorithm must be char/string.');
            end
            if ~isempty(obj.task) && ~any(strcmpi(obj.task, {'classification','regression'}))
                error('predictive_model:InvalidProperty', ...
                    'task must be ''classification'' or ''regression''.');
            end
            if ~isempty(obj.modeloptions)
                if ~iscell(obj.modeloptions)
                    error('predictive_model:InvalidProperty', ...
                        'modeloptions must be a cell array of name/value pairs.');
                end
                % Fixed predicate: option NAMES (odd positions) must be char/string;
                % VALUES (even positions) may be anything.
                names = obj.modeloptions(1:2:end);
                bad = ~cellfun(@(c) ischar(c) || isstring(c), names);
                if any(bad)
                    error('predictive_model:InvalidProperty', ...
                        ['modeloptions: option names (odd positions) must be ' ...
                         'char/string; got non-string in position(s) %s.'], ...
                        num2str(find(bad)*2 - 1));
                end
            end
            if ~isempty(obj.class_labels)
                if ~(iscell(obj.class_labels) && ...
                     all(cellfun(@(c) ischar(c) || isstring(c), obj.class_labels)))
                    error('predictive_model:InvalidProperty', ...
                        'class_labels must be a cell array of strings.');
                end
            end
            if ~isempty(obj.Y_name) && ~(ischar(obj.Y_name) || isstring(obj.Y_name))
                error('predictive_model:InvalidProperty', 'Y_name must be char/string.');
            end
            if ~isempty(obj.X_name) && ~(ischar(obj.X_name) || isstring(obj.X_name))
                error('predictive_model:InvalidProperty', 'X_name must be char/string.');
            end
            if ~isempty(obj.accfun) && ~isa(obj.accfun, 'function_handle')
                error('predictive_model:InvalidProperty', ...
                    'accfun must be a function handle.');
            end

            % Categorised fitted-state structs: must be structs when set.
            cats = {'inputs','cv_partition','fitted_values','weights', ...
                    'weight_stats','error_metrics','descriptions', ...
                    'bootstrap_results','permutation_results', ...
                    'diagnostics','cross_classify','legacy_extras'};
            for i = 1:numel(cats)
                if ~isstruct(obj.(cats{i}))
                    error('predictive_model:InvalidProperty', ...
                        '%s must be a struct.', cats{i});
                end
            end

            % Inputs.Y / inputs.id must be numeric vectors when present.
            if isfield(obj.inputs, 'Y') && ~isempty(obj.inputs.Y)
                validateattributes(obj.inputs.Y, {'numeric'}, {'vector'}, ...
                    mfilename, 'inputs.Y');
            end
            if isfield(obj.inputs, 'id') && ~isempty(obj.inputs.id)
                validateattributes(obj.inputs.id, {'numeric'}, {'vector'}, ...
                    mfilename, 'inputs.id');
            end

            if doverbose
                disp('predictive_model object validated successfully.');
            end
        end


        % -----------------------------------------------------------------
        % LEGACY DEPENDENT ALIASES (read-only)
        % -----------------------------------------------------------------
        % inputs
        function v = get.Y(obj),                            v = predictive_model.field_or_empty(obj.inputs, 'Y'); end
        function v = get.id(obj),                           v = predictive_model.field_or_empty(obj.inputs, 'id'); end
        function v = get.Y_orig(obj),                       v = predictive_model.field_or_empty(obj.inputs, 'Y_orig'); end
        function v = get.INPUTS(obj),                       v = predictive_model.field_or_empty(obj.inputs, 'INPUTS'); end

        % cv_partition
        function v = get.trIdx(obj),                        v = predictive_model.field_or_empty(obj.cv_partition, 'trIdx'); end
        function v = get.teIdx(obj),                        v = predictive_model.field_or_empty(obj.cv_partition, 'teIdx'); end
        function v = get.nfolds(obj),                       v = predictive_model.field_or_empty(obj.cv_partition, 'nfolds'); end
        function v = get.fold_modeloptions(obj),            v = predictive_model.field_or_empty(obj.cv_partition, 'fold_modeloptions'); end
        function v = get.hyperparams_by_fold(obj),          v = predictive_model.field_or_empty(obj.cv_partition, 'hyperparams_by_fold'); end

        % fitted_values
        function v = get.yfit(obj),                         v = predictive_model.field_or_empty(obj.fitted_values, 'yfit'); end
        function v = get.dist_from_hyperplane_xval(obj),    v = predictive_model.field_or_empty(obj.fitted_values, 'dist_from_hyperplane_xval'); end
        function v = get.class_probability_xval(obj),       v = predictive_model.field_or_empty(obj.fitted_values, 'class_probability_xval'); end
        function v = get.scorediff(obj),                    v = predictive_model.field_or_empty(obj.fitted_values, 'scorediff'); end
        function v = get.scores_within_id(obj),             v = predictive_model.field_or_empty(obj.fitted_values, 'scores_within_id'); end
        function v = get.Y_within_id(obj),                  v = predictive_model.field_or_empty(obj.fitted_values, 'Y_within_id'); end
        function v = get.high_vs_low_scores_within_id(obj), v = predictive_model.field_or_empty(obj.fitted_values, 'high_vs_low_scores_within_id'); end
        function v = get.subjfit(obj),                      v = predictive_model.field_or_empty(obj.fitted_values, 'subjfit'); end
        function v = get.predictions(obj),                  v = predictive_model.field_or_empty(obj.fitted_values, 'predictions'); end
        function v = get.trueLabels(obj),                   v = predictive_model.field_or_empty(obj.fitted_values, 'trueLabels'); end

        % weights
        function v = get.w(obj),                            v = predictive_model.field_or_empty(obj.weights, 'w'); end
        function v = get.weight_obj(obj),                   v = predictive_model.field_or_empty(obj.weights, 'weight_obj'); end
        function v = get.mean_vox_weights(obj),             v = predictive_model.field_or_empty(obj.weights, 'mean_vox_weights'); end
        function v = get.vox_weights(obj),                  v = predictive_model.field_or_empty(obj.weights, 'vox_weights'); end
        function v = get.VOXWEIGHTS(obj),                   v = predictive_model.field_or_empty(obj.weights, 'VOXWEIGHTS'); end
        function v = get.my_intercepts(obj),                v = predictive_model.field_or_empty(obj.weights, 'my_intercepts'); end
        function v = get.subjbetas(obj),                    v = predictive_model.field_or_empty(obj.weights, 'subjbetas'); end

        % weight_stats
        function v = get.boot_w(obj),                       v = predictive_model.field_or_empty(obj.weight_stats, 'boot_w'); end
        function v = get.boot_w_ste(obj),                   v = predictive_model.field_or_empty(obj.weight_stats, 'boot_w_ste'); end
        function v = get.boot_w_mean(obj),                  v = predictive_model.field_or_empty(obj.weight_stats, 'boot_w_mean'); end
        function v = get.wZ(obj),                           v = predictive_model.field_or_empty(obj.weight_stats, 'wZ'); end
        function v = get.wP(obj),                           v = predictive_model.field_or_empty(obj.weight_stats, 'wP'); end
        function v = get.wP_fdr_thr(obj),                   v = predictive_model.field_or_empty(obj.weight_stats, 'wP_fdr_thr'); end
        function v = get.boot_w_fdrsig(obj),                v = predictive_model.field_or_empty(obj.weight_stats, 'boot_w_fdrsig'); end
        function v = get.w_thresh_fdr(obj),                 v = predictive_model.field_or_empty(obj.weight_stats, 'w_thresh_fdr'); end

        % error_metrics
        function v = get.crossval_accuracy(obj),                v = predictive_model.field_or_empty(obj.error_metrics, 'crossval_accuracy'); end
        function v = get.crossval_accuracy_within(obj),         v = predictive_model.field_or_empty(obj.error_metrics, 'crossval_accuracy_within'); end
        function v = get.d_singleinterval(obj),                 v = predictive_model.field_or_empty(obj.error_metrics, 'd_singleinterval'); end
        function v = get.d_within(obj),                         v = predictive_model.field_or_empty(obj.error_metrics, 'd_within'); end
        function v = get.classification_d_singleinterval(obj),  v = predictive_model.field_or_empty(obj.error_metrics, 'd_singleinterval'); end
        function v = get.classification_d_within(obj),          v = predictive_model.field_or_empty(obj.error_metrics, 'd_within'); end
        function v = get.regression_d_singleinterval(obj),      v = predictive_model.field_or_empty(obj.error_metrics, 'd_singleinterval'); end
        function v = get.regression_d_within(obj),              v = predictive_model.field_or_empty(obj.error_metrics, 'd_within'); end
        function v = get.prediction_outcome_r(obj),             v = predictive_model.field_or_empty(obj.error_metrics, 'prediction_outcome_r'); end
        function v = get.pred_outcome_r(obj),                   v = predictive_model.field_or_empty(obj.error_metrics, 'prediction_outcome_r'); end
        function v = get.cverr(obj),                            v = predictive_model.field_or_empty(obj.error_metrics, 'cverr'); end
        function v = get.mse(obj),                              v = predictive_model.field_or_empty(obj.error_metrics, 'mse'); end
        function v = get.rmse(obj),                             v = predictive_model.field_or_empty(obj.error_metrics, 'rmse'); end
        function v = get.meanabserr(obj),                       v = predictive_model.field_or_empty(obj.error_metrics, 'meanabserr'); end
        function v = get.r_squared(obj),                        v = predictive_model.field_or_empty(obj.error_metrics, 'r_squared'); end
        function v = get.var_full(obj),                         v = predictive_model.field_or_empty(obj.error_metrics, 'var_full'); end
        function v = get.var_null(obj),                         v = predictive_model.field_or_empty(obj.error_metrics, 'var_null'); end
        function v = get.var_reduction(obj),                    v = predictive_model.field_or_empty(obj.error_metrics, 'var_reduction'); end
        function v = get.pred_err(obj),                         v = predictive_model.field_or_empty(obj.error_metrics, 'pred_err'); end
        function v = get.pred_err_null(obj),                    v = predictive_model.field_or_empty(obj.error_metrics, 'pred_err_null'); end
        function v = get.devs_from_full_model(obj),             v = predictive_model.field_or_empty(obj.error_metrics, 'devs_from_full_model'); end
        function v = get.devs_from_mean_only_model(obj),        v = predictive_model.field_or_empty(obj.error_metrics, 'devs_from_mean_only_model'); end
        function v = get.r_each_subject(obj),                   v = predictive_model.field_or_empty(obj.error_metrics, 'r_each_subject'); end
        function v = get.accuracy(obj),                         v = predictive_model.field_or_empty(obj.error_metrics, 'accuracy'); end
        function v = get.overallAccuracy(obj),                  v = predictive_model.field_or_empty(obj.error_metrics, 'overallAccuracy'); end
        function v = get.err(obj),                              v = predictive_model.field_or_empty(obj.error_metrics, 'err'); end
        function v = get.phi(obj),                              v = predictive_model.field_or_empty(obj.error_metrics, 'phi'); end

        % descriptions
        function v = get.accfun_descrip(obj),                   v = predictive_model.field_or_empty(obj.descriptions, 'accfun_descrip'); end
        function v = get.cverrfun(obj),                         v = predictive_model.field_or_empty(obj.descriptions, 'cverrfun'); end
        function v = get.crossval_accuracy_descrip(obj),        v = predictive_model.field_or_empty(obj.descriptions, 'crossval_accuracy_descrip'); end
        function v = get.prediction_outcome_r_descrip(obj),     v = predictive_model.field_or_empty(obj.descriptions, 'prediction_outcome_r_descrip'); end
        function v = get.pred_err_descrip(obj),                 v = predictive_model.field_or_empty(obj.descriptions, 'pred_err_descrip'); end
        function v = get.var_reduction_descrip(obj),            v = predictive_model.field_or_empty(obj.descriptions, 'var_reduction_descrip'); end
        function v = get.classification_d_within_descrip(obj),  v = predictive_model.field_or_empty(obj.descriptions, 'classification_d_within_descrip'); end
        function v = get.note(obj),                             v = predictive_model.field_or_empty(obj.descriptions, 'note'); end
        function v = get.r_each_subject_note(obj),              v = predictive_model.field_or_empty(obj.descriptions, 'r_each_subject_note'); end
        function v = get.r_each_subject_note2(obj),             v = predictive_model.field_or_empty(obj.descriptions, 'r_each_subject_note2'); end
        function v = get.cvoutput_descrip(obj),                 v = predictive_model.field_or_empty(obj.descriptions, 'cvoutput_descrip'); end
        function v = get.error_type(obj),                       v = predictive_model.field_or_empty(obj.descriptions, 'error_type'); end
        function v = get.function_call(obj),                    v = predictive_model.field_or_empty(obj.descriptions, 'function_call'); end
        function v = get.function_handle(obj),                  v = predictive_model.field_or_empty(obj.descriptions, 'function_handle'); end
        function v = get.algorithm_name(obj),                   v = predictive_model.field_or_empty(obj.descriptions, 'algorithm_name'); end

        % ml_model legacy aliases
        function v = get.ClassificationModel(obj),              v = obj.ml_model; end
        function v = get.SVMModel(obj),                         v = obj.ml_model; end
        function v = get.SVRModel(obj),                         v = obj.ml_model; end

        % bootstrap_results
        function v = get.WTS(obj),                              v = predictive_model.field_or_empty(obj.bootstrap_results, 'WTS'); end

        % diagnostics
        function v = get.mult_obs_within_person(obj),           v = predictive_model.field_or_empty(obj.diagnostics, 'mult_obs_within_person'); end
        function v = get.ROC_forced_choice(obj),                v = predictive_model.field_or_empty(obj.diagnostics, 'ROC_forced_choice'); end
        function v = get.ROC_single_interval(obj),              v = predictive_model.field_or_empty(obj.diagnostics, 'ROC_single_interval'); end
        function v = get.all_reg_hyperparams(obj),              v = predictive_model.field_or_empty(obj.diagnostics, 'all_reg_hyperparams'); end

        % legacy_extras composite fields (brain wrappers)
        function v = get.full_model(obj),                       v = predictive_model.field_or_empty(obj.legacy_extras, 'full_model'); end
        function v = get.data(obj),                             v = predictive_model.field_or_empty(obj.legacy_extras, 'data'); end
        function v = get.covs(obj),                             v = predictive_model.field_or_empty(obj.legacy_extras, 'covs'); end
        function v = get.inputOptions(obj),                     v = predictive_model.field_or_empty(obj.legacy_extras, 'inputOptions'); end

    end


    methods (Access = protected)

        % -----------------------------------------------------------------
        function obj = populate_from_struct(obj, S, doverbose)
            % Walk the input struct and route each field to its categorised
            % home using the static routing table. Unrouted fields land in
            % legacy_extras so nothing is lost.

            routing  = predictive_model.field_routing();
            hp_props = predictive_model.hyperparameter_names();

            % Build a lookup: input-field-name -> {category, subfield}
            keys = routing(:,1);

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

                % 2) top-level fitted-state field (accfun, history, removed_voxels)
                if any(strcmp(name, {'accfun','history','removed_voxels'}))
                    obj.(name) = val;
                    continue
                end

                % 3) categorised routing
                idx = find(strcmp(keys, name), 1);
                if ~isempty(idx)
                    cat = routing{idx, 2};
                    sub = routing{idx, 3};
                    if isempty(sub)
                        % store directly at the category (e.g., ml_model,
                        % fold_models, bootstrap_results when struct, etc.)
                        obj.(cat) = val;
                    else
                        % ensure the category is a struct, then write subfield
                        if ~isstruct(obj.(cat))
                            obj.(cat) = struct();
                        end
                        obj.(cat).(sub) = val;
                    end
                    continue
                end

                % 4) catch-all
                obj.legacy_extras.(name) = val;
                unrouted_to_extras{end+1, 1} = name; %#ok<AGROW>
            end

            % Validate after population.
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
            % Safely read a sub-struct field; return [] if absent.
            if isstruct(s) && isfield(s, name)
                v = s.(name);
            else
                v = [];
            end
        end


        % -----------------------------------------------------------------
        function names = hyperparameter_names()
            % Names of the top-level hyperparameter properties.
            names = {'algorithm','task','modeloptions','random_state', ...
                     'standardize','use_parallel','cv','scorer','nboot', ...
                     'nperm','do_calibrate','Y_name','X_name','class_labels'};
        end


        % -----------------------------------------------------------------
        function routing = field_routing()
            % field_routing  Lookup table from input-struct field name to
            % {category, subfield}. An empty subfield means "store at the
            % category directly".
            %
            % Adding a new wrapper: drop new rows here so the constructor
            % routes their fields without code change elsewhere.

            routing = { ...
                % --- inputs ---
                'Y',                                 'inputs',              'Y'
                'y',                                 'inputs',              'Y'
                'id',                                'inputs',              'id'
                'Y_orig',                            'inputs',              'Y_orig'
                'INPUTS',                            'inputs',              'INPUTS'

                % --- cv_partition ---
                'trIdx',                             'cv_partition',        'trIdx'
                'teIdx',                             'cv_partition',        'teIdx'
                'nfolds',                            'cv_partition',        'nfolds'
                'fold_modeloptions',                 'cv_partition',        'fold_modeloptions'
                'hyperparams_by_fold',               'cv_partition',        'hyperparams_by_fold'
                'cvpartition',                       'cv_partition',        'cvpartition'

                % --- fitted_values ---
                'yfit',                              'fitted_values',       'yfit'
                'dist_from_hyperplane_xval',         'fitted_values',       'dist_from_hyperplane_xval'
                'class_probability_xval',            'fitted_values',       'class_probability_xval'
                'scores_within_id',                  'fitted_values',       'scores_within_id'
                'Y_within_id',                       'fitted_values',       'Y_within_id'
                'scorediff',                         'fitted_values',       'scorediff'
                'high_vs_low_scores_within_id',      'fitted_values',       'high_vs_low_scores_within_id'
                'subjfit',                           'fitted_values',       'subjfit'
                'predictions',                       'fitted_values',       'predictions'
                'trueLabels',                        'fitted_values',       'trueLabels'

                % --- weights ---
                'w',                                 'weights',             'w'
                'mean_vox_weights',                  'weights',             'mean_vox_weights'
                'vox_weights',                       'weights',             'vox_weights'
                'VOXWEIGHTS',                        'weights',             'VOXWEIGHTS'
                'my_intercepts',                     'weights',             'my_intercepts'
                'subjbetas',                         'weights',             'subjbetas'
                'weight_obj',                        'weights',             'weight_obj'

                % --- weight_stats ---
                'boot_w',                            'weight_stats',        'boot_w'
                'boot_w_ste',                        'weight_stats',        'boot_w_ste'
                'boot_w_mean',                       'weight_stats',        'boot_w_mean'
                'wZ',                                'weight_stats',        'wZ'
                'wP',                                'weight_stats',        'wP'
                'wP_fdr_thr',                        'weight_stats',        'wP_fdr_thr'
                'boot_w_fdrsig',                     'weight_stats',        'boot_w_fdrsig'
                'w_thresh_fdr',                      'weight_stats',        'w_thresh_fdr'

                % --- error_metrics (incl. legacy d-prefix translations) ---
                'crossval_accuracy',                 'error_metrics',       'crossval_accuracy'
                'crossval_accuracy_within',          'error_metrics',       'crossval_accuracy_within'
                'd_singleinterval',                  'error_metrics',       'd_singleinterval'
                'd_within',                          'error_metrics',       'd_within'
                'classification_d_singleinterval',   'error_metrics',       'd_singleinterval'
                'regression_d_singleinterval',       'error_metrics',       'd_singleinterval'
                'classification_d_within',           'error_metrics',       'd_within'
                'regression_d_within',               'error_metrics',       'd_within'
                'prediction_outcome_r',              'error_metrics',       'prediction_outcome_r'
                'pred_outcome_r',                    'error_metrics',       'prediction_outcome_r'
                'cverr',                             'error_metrics',       'cverr'
                'mse',                               'error_metrics',       'mse'
                'rmse',                              'error_metrics',       'rmse'
                'meanabserr',                        'error_metrics',       'meanabserr'
                'r_squared',                         'error_metrics',       'r_squared'
                'var_full',                          'error_metrics',       'var_full'
                'var_null',                          'error_metrics',       'var_null'
                'var_reduction',                     'error_metrics',       'var_reduction'
                'pred_err',                          'error_metrics',       'pred_err'
                'pred_err_null',                     'error_metrics',       'pred_err_null'
                'devs_from_full_model',              'error_metrics',       'devs_from_full_model'
                'devs_from_mean_only_model',         'error_metrics',       'devs_from_mean_only_model'
                'r_each_subject',                    'error_metrics',       'r_each_subject'
                'accuracy',                          'error_metrics',       'accuracy'
                'overallAccuracy',                   'error_metrics',       'overallAccuracy'
                'err',                               'error_metrics',       'err'
                'phi',                               'error_metrics',       'phi'

                % --- descriptions ---
                'accfun_descrip',                    'descriptions',        'accfun_descrip'
                'cverrfun',                          'descriptions',        'cverrfun'
                'crossval_accuracy_descrip',         'descriptions',        'crossval_accuracy_descrip'
                'prediction_outcome_r_descrip',      'descriptions',        'prediction_outcome_r_descrip'
                'pred_err_descrip',                  'descriptions',        'pred_err_descrip'
                'var_reduction_descrip',             'descriptions',        'var_reduction_descrip'
                'classification_d_within_descrip',   'descriptions',        'classification_d_within_descrip'
                'note',                              'descriptions',        'note'
                'r_each_subject_note',               'descriptions',        'r_each_subject_note'
                'r_each_subject_note2',              'descriptions',        'r_each_subject_note2'
                'cvoutput_descrip',                  'descriptions',        'cvoutput_descrip'
                'error_type',                        'descriptions',        'error_type'
                'function_call',                     'descriptions',        'function_call'
                'function_handle',                   'descriptions',        'function_handle'
                'algorithm_name',                    'descriptions',        'algorithm_name'

                % --- ml_model (legacy single-object aliases) ---
                'ClassificationModel',               'ml_model',            ''
                'SVMModel',                          'ml_model',            ''
                'SVRModel',                          'ml_model',            ''

                % --- fold_models ---
                'models',                            'fold_models',         ''

                % --- bootstrap_results ---
                'bootstrap',                         'bootstrap_results',   ''
                'WTS',                               'bootstrap_results',   'WTS'
                'boot_weights',                      'bootstrap_results',   'boot_weights'

                % --- diagnostics ---
                'ROC_forced_choice',                 'diagnostics',         'ROC_forced_choice'
                'ROC_single_interval',               'diagnostics',         'ROC_single_interval'
                'mult_obs_within_person',            'diagnostics',         'mult_obs_within_person'
                'all_reg_hyperparams',               'diagnostics',         'all_reg_hyperparams'

                % --- cross_classify ---
                'stats1',                            'cross_classify',      'stats1'
                'stats2',                            'cross_classify',      'stats2'
                'test_results',                      'cross_classify',      'test_results'
                'crosstestfit',                      'cross_classify',      'crosstestfit'
                'cvoutput',                          'cross_classify',      'cvoutput'
                'test_Y',                            'cross_classify',      'test_Y'
                'all',                               'cross_classify',      'all'

                % --- legacy_extras (explicit routing for known wrapper-specific
                %     fields; reads from legacy_extras.<name>) ---
                'covs',                              'legacy_extras',       'covs'
                'data',                              'legacy_extras',       'data'
                'full_model',                        'legacy_extras',       'full_model'
                'inputOptions',                      'legacy_extras',       'inputOptions'
            };
        end

    end

end
