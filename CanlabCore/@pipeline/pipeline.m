classdef pipeline
    % pipeline  Sequential transform steps + a final predictive_model estimator.
    %
    % A scikit-learn-style composition: zero or more *transform steps*
    % (each with a fit_transform / transform contract) feeding a final
    % *estimator* (a @predictive_model, with fit / predict). The whole
    % pipeline behaves like an estimator — fit / predict / score / crossval —
    % so it slots into the same cross-validation machinery.
    %
    % The point of the pipeline is leakage-free cross-validation: when you
    % crossval a pipeline, EVERY step is refit on the training rows of each
    % fold (the held-out rows never touch the transformer fit). This is the
    % correct way to cross-validate PCA-then-regression ("PCR"),
    % standardize-then-SVM, feature-select-then-classify, etc.
    %
    % BUILT-IN TRANSFORM STEPS (referenced by name)
    %   'center'  subtract the training-set feature means
    %   'zscore'  center and divide by training-set feature SDs
    %   'pca'     PCA on the training rows; project onto the top-k PCs.
    %             Param: 'k' (number of components; default = min(n-1, p)).
    %
    % CUSTOM TRANSFORM STEPS
    %   Pass a struct with two function-handle fields:
    %       step.fit_transform = @(Xtr, Ytr) deal(Xtr2, state)
    %       step.transform     = @(state, X)  X2
    %   and, optionally, for weight back-projection to input space:
    %       step.inverse_weights = @(state, w) w_in_input_space
    %
    % :Usage:
    % ::
    %     % PCR: PCA(20) -> ridge regression
    %     est  = predictive_model('algorithm','ridge','task','regression');
    %     pipe = pipeline({ {'pca','k',20} }, est);
    %     pipe = crossval(pipe, X, Y);
    %     pipe.error_metrics.r2.value
    %
    %     % standardize -> linear SVM (classification)
    %     est  = predictive_model('algorithm','svm','task','classification');
    %     pipe = pipeline({'zscore'}, est);
    %     pipe = crossval(pipe, X, Y, 'cv', cv_splitter.stratified_kfold(5));
    %
    % :Inputs (constructor):
    %
    %   **steps:**
    %        a cell array; each element is either a built-in step name
    %        ('center'|'zscore'|'pca'), a cell {name, 'param', value, ...},
    %        or a custom-step struct (see above). {} for no transform.
    %
    %   **estimator:**
    %        a @predictive_model (the final fit/predict stage).
    %
    % :See also:
    %   predictive_model, crossval, cv_splitter, cv_scorer
    %
    % -------------------------------------------------------------------------
    %     Copyright (C) 2025  CANlab — GNU GPL v3 or later.
    % -------------------------------------------------------------------------

    properties
        steps        = {};      % cell of normalized step structs
        estimator    = [];      % final @predictive_model
        cv           = [];      % cv_splitter (optional; defaults at crossval)
        scorer       = [];      % cv_scorer  (optional)
        random_state = [];
    end

    properties (SetAccess = protected)
        fit_type      = '';
        fitted_values = struct();
        error_metrics = struct();
        cv_partition  = struct();
        weights       = struct();   % .w in input space (back-projected)
        history       = {};
    end

    methods

        % -----------------------------------------------------------------
        function obj = pipeline(steps, estimator)
            if nargin == 0, return; end
            if nargin >= 1 && ~isempty(steps)
                if ~iscell(steps), steps = {steps}; end
                obj.steps = cellfun(@pipeline.normalize_step, steps, 'uniform', false);
            end
            if nargin >= 2
                if ~isa(estimator, 'predictive_model')
                    error('pipeline:BadEstimator', ...
                        'estimator must be a @predictive_model.');
                end
                obj.estimator = estimator;
            end
        end


        % -----------------------------------------------------------------
        function obj = fit(obj, X, Y, varargin)
            % fit  Fit every transform step then the estimator (in-sample).
            if isempty(obj.estimator)
                error('pipeline:NoEstimator', 'Set obj.estimator before fitting.');
            end
            if ~isempty(obj.random_state), rng(obj.random_state); end

            Xc = X;
            for i = 1:numel(obj.steps)
                [Xc, obj.steps{i}] = pipeline.step_fit_transform(obj.steps{i}, Xc, Y);
            end
            obj.estimator = fit(obj.estimator, Xc, Y, varargin{:});

            obj.weights = obj.back_project_weights();
            obj.fit_type = 'insample';
            obj.history{end+1, 1} = sprintf('fit: %d step(s) + %s on %d obs', ...
                numel(obj.steps), char(string(obj.estimator.algorithm)), size(X, 1));
        end


        % -----------------------------------------------------------------
        function [yhat, score_out] = predict(obj, X)
            % predict  Transform through every step, then estimator.predict.
            Xc = X;
            for i = 1:numel(obj.steps)
                Xc = pipeline.step_transform(obj.steps{i}, Xc);
            end
            if nargout > 1
                [yhat, score_out] = predict(obj.estimator, Xc);
            else
                yhat = predict(obj.estimator, Xc);
            end
        end


        % -----------------------------------------------------------------
        function v = score(obj, X, Y)
            % score  Evaluate obj.scorer (or estimator default) on predictions.
            sc = obj.scorer;
            if isempty(sc), sc = obj.estimator.scorer; end
            if isempty(sc), sc = cv_scorer.default_for_task(obj.estimator.task); end
            if sc.needs_continuous
                [yhat, scores] = predict(obj, X);
                v = sc.score(Y, yhat, scores);
            else
                yhat = predict(obj, X);
                v = sc.score(Y, yhat);
            end
        end


        % -----------------------------------------------------------------
        function new_obj = clone(obj)
            % clone  Fresh, unfitted copy: steps stripped of learned state,
            % estimator cloned. Hyperparameters preserved.
            new_obj = pipeline();
            new_obj.estimator    = clone(obj.estimator);
            new_obj.cv           = obj.cv;
            new_obj.scorer       = obj.scorer;
            new_obj.random_state = obj.random_state;
            new_obj.steps        = obj.steps;
            for i = 1:numel(new_obj.steps)
                new_obj.steps{i}.state = [];
            end
        end


        % -----------------------------------------------------------------
        function tf = is_fitted(obj)
            tf = ~isempty(obj.fit_type) || ...
                (~isempty(obj.estimator) && is_fitted(obj.estimator));
        end


        % -----------------------------------------------------------------
        function w = back_project_weights(obj)
            % back_project_weights  Map the estimator's coefficient vector
            % from the transformed feature space back to the pipeline input
            % space by walking the steps in reverse. Returns struct .w (or
            % empty .w if any step cannot be inverted / estimator is
            % nonlinear).
            w = struct('w', []);
            if isempty(obj.estimator) || ~isfield(obj.estimator.weights, 'w') ...
                    || isempty(obj.estimator.weights.w)
                return
            end
            wv = obj.estimator.weights.w(:);

            % Expand to the FULL transformed-feature length if the estimator
            % dropped any transformed features (omitted_features lives in
            % the estimator's transformed input space, not voxel space).
            of = obj.estimator.omitted_features;
            if islogical(of) && any(of) && numel(of) == numel(wv) + sum(of)
                full = zeros(numel(of), 1);
                full(~of) = wv;
                wv = full;
            end

            for i = numel(obj.steps):-1:1
                wv = pipeline.step_inverse_weights(obj.steps{i}, wv);
                if isempty(wv), w = struct('w', []); return; end
            end
            w.w = wv;
        end


        % -----------------------------------------------------------------
        function si = weight_image(obj, source, varargin)
            % weight_image  Back-project weights to voxel space as a
            % @statistic_image (so montage(weight_image(pipe,source)) works,
            % even for PCA-reduced PCR pipelines). Delegates to the
            % estimator's weight_image after stashing the back-projected
            % weights into a temporary predictive_model.
            wstruct = obj.back_project_weights();
            if isempty(wstruct.w)
                error('pipeline:weight_image:NoInputWeights', ...
                    ['Could not back-project weights to input space (a step ' ...
                     'is non-invertible or the estimator is nonlinear).']);
            end
            % Back-projected weights are already in the full input (voxel)
            % space, so weight_image direct-matches them to `source`; pass
            % empty omitted_features.
            tmp = pipeline.inject_weights(clone(obj.estimator), wstruct.w, []);
            si = weight_image(tmp, source, varargin{:});
        end

    end


    methods (Static)

        % -----------------------------------------------------------------
        function s = normalize_step(stepspec)
            % normalize_step  Coerce a step spec into a canonical struct.
            s = struct('name', '', 'params', struct(), 'state', [], ...
                       'is_custom', false, 'fit_transform', [], ...
                       'transform', [], 'inverse_weights', []);

            if ischar(stepspec) || isstring(stepspec)
                s.name = char(stepspec);

            elseif iscell(stepspec)
                s.name = char(stepspec{1});
                nv = stepspec(2:end);
                for k = 1:2:numel(nv)
                    s.params.(nv{k}) = nv{k+1};
                end

            elseif isstruct(stepspec)
                if ~isfield(stepspec, 'fit_transform') || ~isfield(stepspec, 'transform')
                    error('pipeline:BadStep', ...
                        'Custom step structs must have fit_transform and transform fields.');
                end
                s.is_custom     = true;
                s.name          = pipeline.field_or(stepspec, 'name', 'custom');
                s.fit_transform = stepspec.fit_transform;
                s.transform     = stepspec.transform;
                s.inverse_weights = pipeline.field_or(stepspec, 'inverse_weights', []);

            else
                error('pipeline:BadStep', ...
                    'Step must be a name, a {name, nv...} cell, or a custom struct.');
            end

            if ~s.is_custom && ~any(strcmpi(s.name, {'center','zscore','pca'}))
                error('pipeline:UnknownStep', ...
                    'Unknown built-in step ''%s''. Use center | zscore | pca, or a custom struct.', s.name);
            end
        end


        % -----------------------------------------------------------------
        function [X2, s] = step_fit_transform(s, Xtr, Ytr)
            if s.is_custom
                [X2, s.state] = s.fit_transform(Xtr, Ytr);
                return
            end
            switch lower(s.name)
                case 'center'
                    s.state.mu = mean(Xtr, 1);
                    X2 = Xtr - s.state.mu;
                case 'zscore'
                    s.state.mu = mean(Xtr, 1);
                    sd = std(Xtr, 0, 1); sd(sd == 0) = 1;
                    s.state.sd = sd;
                    X2 = (Xtr - s.state.mu) ./ sd;
                case 'pca'
                    s.state.mu = mean(Xtr, 1);
                    Xc = Xtr - s.state.mu;
                    [coeff, scoreM] = pca_econ(Xc);
                    k = pipeline.field_or(s.params, 'k', size(coeff, 2));
                    k = max(1, min(k, size(coeff, 2)));
                    s.state.coeff = coeff(:, 1:k);
                    X2 = scoreM(:, 1:k);
            end
        end


        % -----------------------------------------------------------------
        function X2 = step_transform(s, X)
            if s.is_custom
                X2 = s.transform(s.state, X);
                return
            end
            switch lower(s.name)
                case 'center'
                    X2 = X - s.state.mu;
                case 'zscore'
                    X2 = (X - s.state.mu) ./ s.state.sd;
                case 'pca'
                    X2 = (X - s.state.mu) * s.state.coeff;
            end
        end


        % -----------------------------------------------------------------
        function w_in = step_inverse_weights(s, w)
            % Map weights from this step's OUTPUT space back to its INPUT
            % space. Returns [] if the step is non-invertible.
            if s.is_custom
                if isempty(s.inverse_weights)
                    w_in = [];
                else
                    w_in = s.inverse_weights(s.state, w);
                end
                return
            end
            switch lower(s.name)
                case 'center'
                    w_in = w;                         % centering: weights unchanged
                case 'zscore'
                    w_in = w ./ s.state.sd(:);        % undo the per-feature scaling
                case 'pca'
                    w_in = s.state.coeff * w;         % coeff: [p x k], w: [k x 1]
                otherwise
                    w_in = [];
            end
        end


        % -----------------------------------------------------------------
        function v = field_or(s, name, default)
            if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
                v = s.(name);
            else
                v = default;
            end
        end


        % -----------------------------------------------------------------
        function pm = inject_weights(pm, w, omitted_features)
            % Stash a weight vector into a predictive_model so weight_image
            % can map it (used by pipeline.weight_image). Uses the public
            % struct-constructor route so SetAccess=protected is respected.
            S = struct('w', w(:), 'omitted_features', omitted_features, ...
                       'algorithm', char(string(pm.algorithm)), ...
                       'fit_type', 'insample');
            pm2 = predictive_model(S, 'noverbose');
            pm2.algorithm = pm.algorithm;
            pm2.task      = pm.task;
            pm = pm2;
        end

    end

end


% ----------------------------- local helper -------------------------------
function [coeff, scoreM] = pca_econ(Xc)
% Economy PCA via SVD on already-centered Xc ([n x p]). Returns loadings
% coeff ([p x r]) and scores scoreM ([n x r]), r = rank.
    [U, S, V] = svd(Xc, 'econ');
    coeff  = V;
    scoreM = U * S;
end
