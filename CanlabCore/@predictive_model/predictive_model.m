classdef predictive_model
    % predictive_model Object representing a predictive model.
    %
    % :Usage:
    % ::
    %     pm_obj = predictive_model();       % Create an empty predictive_model object.
    %     pm_obj = predictive_model(pm);     % Create a predictive_model object from structure pm.
    %
    % :Properties:
    %
    %   **Y:**
    %        [1000×1 double] Actual (observed) outcomes (for SVM, 1 or -1).
    %
    %   **id:**
    %        [1000×1 double] Grouping variable for within-participant observations.
    %
    %   **class_labels:**
    %        Cell array with strings specifying class names, e.g., {'No pain' 'Pain'}
    %
    %   **Y_name:**
    %        String specifying name of outcome variable, e.g., {'Pain_at_MRI'}
    %
    %   **X_name:**
    %        String specifying name of predictor variable, e.g., {'fMRI pattern response'}
    %
    %   **modeloptions:**
    %        Cell array with strings specifying model options, e.g., {'KernelFunction','linear'}.
    %
    %   **accfun:**
    %        Function handle to compute accuracy, e.g., @(Y,yfit)100.*nansum(Y==yfit)./sum(~isnan(Y)).
    %
    %   **trIdx:**
    %        {1×10 cell} Cell array of logical vectors for training indices.
    %
    %   **teIdx:**
    %        {1×10 cell} Cell array of logical vectors for testing indices.
    %
    %   **nfolds:**
    %        Scalar (e.g., 10). Number of folds in the holdout set.
    %
    %   **dist_from_hyperplane_xval:**
    %        [1000×1 double] Cross-validated distances from the hyperplane.
    %
    %   **yfit:**
    %        [1000×1 double] Predicted outcomes.
    %
    %   **class_probability_xval:**
    %        [1000×1 double] Cross-validated probabilities for Class 1 (using Platt scaling).
    %
    %   **crossval_accuracy:**
    %        Scalar (e.g., 55.5000). Cross-validated accuracy.
    %
    %   **classification_d_singleinterval:**
    %        Scalar (e.g., 0.1847). Classification effect size.
    %
    %   **mult_obs_within_person:**
    %        Scalar (e.g., 0). Flag for multiple observations within a person.
    %
    %   **ClassificationModel:**
    %        Model object (e.g., a ClassificationSVM) trained on the full dataset.
    %
    %   **w:**
    %        Numeric vector ([148×1 double]) of model weights/betas.
    %
    %   **boot_w:**
    %        Numeric vector of bootstrapped weights.
    %
    %   **boot_w_ste:**
    %        Numeric vector of bootstrapped standard errors for model weights.
    %
    %   **boot_w_mean:**
    %        Numeric vector of bootstrapped mean model weights.
    %
    %   **wZ:**
    %        Numeric vector of bootstrapped Z-scores for model weights.
    %
    %   **wP:**
    %        Numeric vector of bootstrapped P-values for individual model weights.
    %
    %   **wP_fdr_thr:**
    %        Numeric value (scalar or vector) for the FDR threshold (q < 0.05).
    %
    %   **boot_w_fdrsig:**
    %        Logical vector indicating which bootstrapped weights are FDR-significant.
    %
    %   **w_thresh_fdr:**
    %        Numeric vector of thresholded weights at FDR q < 0.05.
    %
    %   **accfun_descrip:**
    %        A string describing the accuracy function (accfun) used by the model.
    %
    %   **cverrfun:**
    %        Function handle for computing cross-validation error.
    %
    %   **cverr:**
    %        Numeric vector containing cross-validation error values.
    %
    %   **crossval_accuracy_descrip:**
    %        A string description of the cross-validated accuracy metric.
    %
    %   **prediction_outcome_r:**
    %        Numeric scalar representing the correlation (r) between predicted and actual outcomes.
    %
    %   **prediction_outcome_r_descrip:**
    %        A string description of the prediction-outcome correlation metric.
    %
    %   **d_singleinterval:**
    %        Numeric scalar representing the effect size for regression analysis on a single interval.
    %
    %   **d_within:**
    %        Numeric scalar representing the within-person classification effect size.
    %
    %
    % :Methods (forthcoming):
    %   validate_object   - Validates the object properties.
    %   train             - Train the predictive model. (forthcoming)
    %   test              - Test the predictive model. (forthcoming)
    %   crossval          - Cross-validation of the model. (forthcoming)
    %   report            - Generate a performance report. (forthcoming)
    %   plot              - Plot predictions versus outcomes. (forthcoming)
    %   montage           - Display a montage of images. (forthcoming)
    %   select_features   - Feature selection. (forthcoming)
    %   bootstrap         - Bootstrap analysis and reproducibility. (forthcoming)
    %   error_analysis    - Analyze misclassified images and errors. (forthcoming)
    %   permutation_test  - Permutation testing for significance. (forthcoming)
    %   confusionchart    - Create a confusion chart. (forthcoming)
    %   rocplot           - Plot ROC curves. (forthcoming)
    %
    % :Examples:
    % ::
    %     load('pm.mat');  % Assume pm is a structure with the required fields.
    %     pm_obj = predictive_model(pm);
    %     pm_obj.validate_object();  % Validate the predictive_model object.
    %
    % :References:
    %   (Add relevant citations here.)
    %
    % :See also:
    %   validate_object
    %
    % -------------------------------------------------------------------------
    %     Author and copyright information:
    %
    %     Copyright (C) 2025  Your Name
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

    properties
        Y                       % [1000×1 double]
        id                      % [1000×1 double]
        class_labels = {}; 
        Y_name = {};
        X_name = {};
        modeloptions            % Cell array, e.g., {'KernelFunction','linear'}
        accfun                  % Function handle, e.g., @(Y,yfit)100.*nansum(Y==yfit)./sum(~isnan(Y))
        trIdx        = {};           % {1×10 cell} of logical vectors
        teIdx        = {};            % {1×10 cell} of logical vectors
        nfolds                  % Scalar, e.g., 10
        dist_from_hyperplane_xval  % [1000×1 double]
        yfit                    % [1000×1 double]
        class_probability_xval  % [1000×1 double]
        crossval_accuracy       % Scalar, e.g., 55.5000
        classification_d_singleinterval  % Scalar, e.g., 0.1847
        mult_obs_within_person  % Scalar, e.g., 0
        ClassificationModel     % Model object, e.g., [1×1 ClassificationSVM]
        w                       % Numeric vector, e.g., [148×1 double]
        boot_w                  % Numeric vector
        boot_w_ste              % Numeric vector
        boot_w_mean             % Numeric vector
        wZ                      % Numeric vector
        wP                      % Numeric vector
        wP_fdr_thr              % Numeric (scalar or vector)
        boot_w_fdrsig           % Logical vector
        w_thresh_fdr            % Numeric vector
        d_singleinterval
        d_within
        accfun_descrip
        cverrfun
        cverr
        crossval_accuracy_descrip
        prediction_outcome_r
        prediction_outcome_r_descrip

    end

    methods
        function obj = predictive_model(varargin)
            % predictive_model Constructs a predictive_model object.
            %
            % :Usage:
            % ::
            %     pm_obj = predictive_model();       % Creates an empty predictive_model object.
            %     pm_obj = predictive_model(pm);     % Creates a predictive_model object from structure pm.
            %     pm_obj = predictive_model(pm, 'noverbose'); % Creates object with verbose output suppressed.
            %
            % If a structure pm is provided, only fields that match the predefined
            % properties are assigned; additional fields are discarded.

            % Default verbose flag.
            doverbose = true;

            % Check for optional 'noverbose' flag.
            if nargin > 1
                for i = 1:length(varargin)
                    if ischar(varargin{i}) && strcmpi(varargin{i}, 'noverbose')
                        doverbose = false;
                        break;
                    end
                end
            end

            % No arguments: return an empty object.
            if nargin == 0
                return;
            end

            % If a structure is provided, assign matching fields.
            if nargin >= 1 && isstruct(varargin{1})
                pm = varargin{1};
                propNames = properties(obj);
                for i = 1:length(propNames)
                    if isfield(pm, propNames{i})
                        obj.(propNames{i}) = pm.(propNames{i});
                    else
                        obj.(propNames{i}) = [];  % Allow properties to be empty
                    end
                end

            else
                error('predictive_model:InvalidInput', 'Input must be a structure.');
            end

            % Special assignments to accommodate legacy xval_SVM and xval_SVR
            % -------------------------------------------------------------
            if isfield(pm, 'SVMModel') && isobject(pm.SVMModel)
                obj.ClassificationModel = pm.SVMModel;
                pm = rmfield(pm, 'SVMModel'); % for reporting purposes, to prevent reporting as non-copied field
            end

            if isfield(pm, 'SVRModel') && isobject(pm.SVRModel)
                obj.ClassificationModel = pm.SVRModel;
                pm = rmfield(pm, 'SVRModel'); % for reporting purposes, to prevent reporting as non-copied field
            end

            if isfield(pm, 'classification_d_singleinterval') && ~isempty(pm.classification_d_singleinterval)
                obj.d_singleinterval = pm.classification_d_singleinterval;
                pm = rmfield(pm, 'classification_d_singleinterval'); % for reporting purposes, to prevent reporting as non-copied field
            end

            if isfield(pm, 'regression_d_singleinterval') && ~isempty(pm.regression_d_singleinterval)
                obj.d_singleinterval = pm.regression_d_singleinterval;
                pm = rmfield(pm, 'regression_d_singleinterval'); % for reporting purposes, to prevent reporting as non-copied field
            end

            if isfield(pm, 'classification_d_within') && ~isempty(pm.classification_d_within)
                obj.d_within = pm.classification_d_within;
                pm = rmfield(pm, 'classification_d_within'); % for reporting purposes, to prevent reporting as non-copied field
            end

            if isfield(pm, 'regression_d_within') && ~isempty(pm.regression_d_within)
                obj.d_within = pm.regression_d_within;
                pm = rmfield(pm, 'regression_d_within'); % for reporting purposes, to prevent reporting as non-copied field
            end


            % Validate the object.
            obj.validate_object(varargin{:});


            % Report non-copied fields, if any

            if nargin >= 1 && isstruct(varargin{1}) && doverbose
                % If extra fields exist, print them.
                fn = fieldnames(pm);
                extraFields = {};

                for i = 1:length(fn)

                    if ~ismember(fn{i}, propNames)
                        extraFields{end + 1, 1} = fn{i};
                    end

                end

                if ~isempty(extraFields)

                    fprintf('The following fields in the input structure were not copied into the object properties:\n');
                    disp(extraFields);

                end
            end % doverbose print extra fields

        end


        function validate_object(obj, varargin)
            % validate_object Validates the data types for predictive_model properties.
            %
            % :Usage:
            % ::
            %     pm_obj.validate_object()
            %
            % This method uses anonymous function handles and validateattributes to check:
            %   - Y and id are numeric vectors.
            %   - modeloptions is a cell array of strings.
            %   - accfun is a function handle.
            %   - trIdx and teIdx are cell arrays.
            %   - nfolds is a numeric scalar.
            %   - dist_from_hyperplane_xval, yfit, class_probability_xval are numeric vectors.
            %   - crossval_accuracy, classification_d_singleinterval, mult_obs_within_person are numeric scalars.
            %   - ClassificationModel is an object.
            %   - w, boot_w, boot_w_ste, boot_w_mean, wZ, wP, w_thresh_fdr are numeric vectors.
            %   - wP_fdr_thr is numeric.
            %   - boot_w_fdrsig is a logical vector.

            doverbose = true;
            if any(strcmp(varargin, 'noverbose')), doverbose = false; end

            % Define anonymous validation functions.
            validateNumericVector = @(x, name) validateattributes(x, {'numeric'}, {'vector'}, mfilename, name);
            validateNumericScalar = @(x, name) validateattributes(x, {'numeric'}, {'scalar'}, mfilename, name);
            validateLogicalScalar = @(x, name) validateattributes(x, {'logical'}, {'scalar'}, mfilename, name);
            % validateCellOfStrings = @(x, name) (iscell(x) && all(cellfun(@(c) ischar(c) || isstring(c), x))) || ...
            %     error('predictive_model:InvalidProperty', '%s must be a cell array of strings.', name);
            % This line expects modeloptions to be a cell array of strings, but it should actually accept key-value pairs (alternating strings and values), 
            % as required by SVM functions. To fix this, modify the validation logic as below to allow cell arrays containing both strings (option names) and 
            % their corresponding values, or temporarily comment out the validation line for modeloptions. Byeol 2025/06/24
            validateModeloptionsCell = @(x, name) (iscell(x) && all(cellfun(@(c,idx) mod(idx,2)==1 && (ischar(c)||isstring(c)) || mod(idx,2)==0, x, num2cell(1:numel(x))))) || ...
                error('predictive_model:InvalidProperty', '%s must be a cell array with alternating strings and values.', name);

            validateFunctionHandle = @(x, name) validateattributes(x, {'function_handle'}, {}, mfilename, name);
            validateCell = @(x, name) validateattributes(x, {'cell'}, {}, mfilename, name);
            validateChar = @(x, name) validateattributes(x, {'char'}, {}, mfilename, name);

            % Validate Y
            if ~isempty(obj.Y)
                validateNumericVector(obj.Y, 'Y');
            end

            % Validate id
            if ~isempty(obj.id)
                validateNumericVector(obj.id, 'id');
            end

            if ~isempty(obj.class_labels)
                validateCellOfStrings(obj.class_labels, 'class_labels');
            end

            if ~isempty(obj.Y_name)
                validateChar(obj.Y_name, 'Y_name');
            end

            if ~isempty(obj.X_name)
                validateChar(obj.X_name, 'X_name');
            end

            % Validate modeloptions
            if ~isempty(obj.modeloptions)
                validateModeloptionsCell(obj.modeloptions, 'modeloptions');
            end

            % Validate accfun
            if ~isempty(obj.accfun)
                validateFunctionHandle(obj.accfun, 'accfun');
            end

            % Validate trIdx and teIdx
            if ~isempty(obj.trIdx)
                validateCell(obj.trIdx, 'trIdx');
            end
            if ~isempty(obj.teIdx)
                validateCell(obj.teIdx, 'teIdx');
            end

            % Validate nfolds
            if ~isempty(obj.nfolds)
                validateNumericScalar(obj.nfolds, 'nfolds');
            end

            % Validate dist_from_hyperplane_xval, yfit, class_probability_xval
            if ~isempty(obj.dist_from_hyperplane_xval)
                validateNumericVector(obj.dist_from_hyperplane_xval, 'dist_from_hyperplane_xval');
            end
            if ~isempty(obj.yfit)
                validateNumericVector(obj.yfit, 'yfit');
            end
            if ~isempty(obj.class_probability_xval)
                validateNumericVector(obj.class_probability_xval, 'class_probability_xval');
            end

            % Validate crossval_accuracy, classification_d_singleinterval, mult_obs_within_person
            if ~isempty(obj.crossval_accuracy)
                validateNumericScalar(obj.crossval_accuracy, 'crossval_accuracy');
            end
            if ~isempty(obj.classification_d_singleinterval)
                validateNumericScalar(obj.classification_d_singleinterval, 'classification_d_singleinterval');
            end
            if ~isempty(obj.mult_obs_within_person)
                validateLogicalScalar(obj.mult_obs_within_person, 'mult_obs_within_person');
            end

            % Validate ClassificationModel (allow empty or any object)
            if ~isempty(obj.ClassificationModel) && ~isobject(obj.ClassificationModel)
                error('predictive_model:InvalidProperty', 'ClassificationModel must be an object.');
            end

            % Validate w, boot_w, boot_w_ste, boot_w_mean, wZ, wP, w_thresh_fdr
            if ~isempty(obj.w)
                validateNumericVector(obj.w, 'w');
            end
            if ~isempty(obj.boot_w)
                validateNumericVector(obj.boot_w, 'boot_w');
            end
            if ~isempty(obj.boot_w_ste)
                validateNumericVector(obj.boot_w_ste, 'boot_w_ste');
            end
            if ~isempty(obj.boot_w_mean)
                validateNumericVector(obj.boot_w_mean, 'boot_w_mean');
            end
            if ~isempty(obj.wZ)
                validateNumericVector(obj.wZ, 'wZ');
            end
            if ~isempty(obj.wP)
                validateNumericVector(obj.wP, 'wP');
            end
            if ~isempty(obj.w_thresh_fdr)
                validateNumericVector(obj.w_thresh_fdr, 'w_thresh_fdr');
            end

            % Validate wP_fdr_thr (allow scalar or vector)
            if ~isempty(obj.wP_fdr_thr)
                validateattributes(obj.wP_fdr_thr, {'numeric'}, {}, mfilename, 'wP_fdr_thr');
            end

            % Validate boot_w_fdrsig
            if ~isempty(obj.boot_w_fdrsig)
                validateattributes(obj.boot_w_fdrsig, {'logical'}, {'vector'}, mfilename, 'boot_w_fdrsig');
            end

            % Validate accfun_descrip (should be a char array or string)
            if ~isempty(obj.accfun_descrip)
                validateattributes(obj.accfun_descrip, {'char','string'}, {}, mfilename, 'accfun_descrip');
            end

            % Validate cverrfun (must be a function handle)
            if ~isempty(obj.cverrfun)
                validateattributes(obj.cverrfun, {'function_handle'}, {}, mfilename, 'cverrfun');
            end

            % Validate cverr (must be a numeric vector)
            if ~isempty(obj.cverr)
                validateattributes(obj.cverr, {'char','string'}, {}, mfilename, 'cverr');
            end

            % Validate crossval_accuracy_descrip (should be a char array or string)
            if ~isempty(obj.crossval_accuracy_descrip)
                validateattributes(obj.crossval_accuracy_descrip, {'char','string'}, {}, mfilename, 'crossval_accuracy_descrip');
            end

            % Validate prediction_outcome_r (should be a numeric scalar)
            if ~isempty(obj.prediction_outcome_r)
                validateattributes(obj.prediction_outcome_r, {'numeric'}, {'scalar'}, mfilename, 'prediction_outcome_r');
            end

            % Validate prediction_outcome_r_descrip (should be a char array or string)
            if ~isempty(obj.prediction_outcome_r_descrip)
                validateattributes(obj.prediction_outcome_r_descrip, {'char','string'}, {}, mfilename, 'prediction_outcome_r_descrip');
            end

            % Validate d_singleinterval (numeric scalar)
            if ~isempty(obj.d_singleinterval)
                validateattributes(obj.d_singleinterval, {'numeric'}, {'scalar'}, mfilename, 'regression_d_singleinterval');
            end

            % Validate d_within (numeric scalar)
            if ~isempty(obj.d_within)
                validateattributes(obj.d_within, {'numeric'}, {'scalar'}, mfilename, 'd_within');
            end

            if doverbose

                disp('predictive_model object validated successfully.');

            end


        end

        % Forthcoming methods:
        function train(obj)
            % train Train the predictive model.
            % (forthcoming)
            error('Method train is forthcoming.');
        end

        function test(obj)
            % test Test the predictive model.
            % (forthcoming)
            error('Method test is forthcoming.');
        end

        function crossval(obj)
            % crossval Perform cross-validation on the predictive model.
            % (forthcoming)
            error('Method crossval is forthcoming.');
        end

        function report(obj)
            % report Generate a performance report for the predictive model.
            % (forthcoming)
            error('Method report is forthcoming.');
        end

        function plot(obj)
            % plot Plot predictions versus outcomes.
            % (forthcoming)
            error('Method plot is forthcoming.');
        end

        function montage(obj)
            % montage Display a montage of relevant images.
            % (forthcoming)
            error('Method montage is forthcoming.');
        end

        function select_features(obj)
            % select_features Perform feature selection.
            % (forthcoming)
            error('Method select_features is forthcoming.');
        end

        function bootstrap(obj)
            % bootstrap Perform bootstrap analysis and reproducibility assessment.
            % (forthcoming)
            error('Method bootstrap is forthcoming.');
        end

        function error_analysis(obj)
            % error_analysis Analyze misclassified images and error metrics.
            % (forthcoming)
            error('Method error_analysis is forthcoming.');
        end

        function permutation_test(obj)
            % permutation_test Perform permutation testing for statistical significance.
            % (forthcoming)
            error('Method permutation_test is forthcoming.');
        end

        function confusionchart(obj)
            % confusionchart Create a confusion chart of model performance.
            % (forthcoming)
            error('Method confusionchart is forthcoming.');
        end

        function rocplot(obj)
            % rocplot Plot ROC curves for the predictive model.
            % (forthcoming)
            error('Method rocplot is forthcoming.');
        end
    end
end