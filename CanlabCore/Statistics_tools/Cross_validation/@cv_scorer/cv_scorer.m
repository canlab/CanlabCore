classdef cv_scorer
    % cv_scorer  Sklearn-style scorer object.
    %
    % Wraps a scoring function with metadata so cross-validation and
    % hyperparameter search can argmax uniformly (the
    % "greater-is-better" invariant). For loss-like metrics (MSE, MAE,
    % log-loss), set greater_is_better = false; cross-validation code
    % negates internally so a single argmax loop works.
    %
    % :Usage:
    % ::
    %    s = cv_scorer.accuracy();
    %    v = s.score(Y, yhat);                  % discrete predictions
    %    s = cv_scorer.roc_auc();
    %    v = s.score(Y, yhat, scores);          % needs continuous scores
    %    s = cv_scorer.make('rmse');            % look up by name
    %    list = cv_scorer.names();              % registered names
    %
    % :Custom scorer:
    % ::
    %    s = cv_scorer('my_metric', @(Y,yhat,~) ..., 'greater_is_better', false);
    %
    % Built-in scorers:
    %   Classification: accuracy, balanced_accuracy, f1, roc_auc, log_loss
    %   Regression:     mse, rmse, mae, r2, pearson_r

    properties
        name             = '';
        fn               = [];      % function handle: @(Y, yhat, scores) -> float
        greater_is_better = true;
        needs_continuous = false;   % true if scorer needs continuous scores (3rd arg)
    end

    methods
        % -----------------------------------------------------------------
        function obj = cv_scorer(name, fn, varargin)
            if nargin == 0, return; end
            obj.name = name;
            if nargin >= 2, obj.fn = fn; end
            p = inputParser;
            addParameter(p, 'greater_is_better', obj.greater_is_better);
            addParameter(p, 'needs_continuous',  obj.needs_continuous);
            parse(p, varargin{:});
            obj.greater_is_better = p.Results.greater_is_better;
            obj.needs_continuous  = p.Results.needs_continuous;
        end

        % -----------------------------------------------------------------
        function v = score(obj, Y, yhat, scores)
            if nargin < 4, scores = []; end
            Y = Y(:); yhat = yhat(:);
            v = obj.fn(Y, yhat, scores);
        end
    end


    methods (Static)
        % -------------------- Factory functions --------------------

        % --- Classification ---
        function s = accuracy()
            s = cv_scorer('accuracy', @(Y,yh,~) mean(Y == yh), 'greater_is_better', true);
        end

        function s = balanced_accuracy()
            s = cv_scorer('balanced_accuracy', @cv_scorer.fn_balanced_accuracy, ...
                          'greater_is_better', true);
        end

        function s = f1()
            s = cv_scorer('f1', @cv_scorer.fn_f1_binary, 'greater_is_better', true);
        end

        function s = roc_auc()
            s = cv_scorer('roc_auc', @cv_scorer.fn_roc_auc, ...
                          'greater_is_better', true, 'needs_continuous', true);
        end

        function s = log_loss()
            s = cv_scorer('log_loss', @cv_scorer.fn_log_loss, ...
                          'greater_is_better', false, 'needs_continuous', true);
        end

        % --- Regression ---
        function s = mse()
            s = cv_scorer('mse', @(Y,yh,~) mean((Y-yh).^2), 'greater_is_better', false);
        end

        function s = rmse()
            s = cv_scorer('rmse', @(Y,yh,~) sqrt(mean((Y-yh).^2)), 'greater_is_better', false);
        end

        function s = mae()
            s = cv_scorer('mae', @(Y,yh,~) mean(abs(Y-yh)), 'greater_is_better', false);
        end

        function s = r2()
            s = cv_scorer('r2', @cv_scorer.fn_r2, 'greater_is_better', true);
        end

        function s = pearson_r()
            s = cv_scorer('pearson_r', @(Y,yh,~) cv_scorer.safe_corr(Y, yh), ...
                          'greater_is_better', true);
        end

        % -------------------- Registry lookup --------------------
        function s = make(name)
            % Look up a scorer by registered name.
            switch lower(name)
                case 'accuracy',          s = cv_scorer.accuracy();
                case 'balanced_accuracy', s = cv_scorer.balanced_accuracy();
                case 'f1',                s = cv_scorer.f1();
                case 'roc_auc',           s = cv_scorer.roc_auc();
                case 'log_loss',          s = cv_scorer.log_loss();
                case 'mse',               s = cv_scorer.mse();
                case 'rmse',              s = cv_scorer.rmse();
                case 'mae',               s = cv_scorer.mae();
                case 'r2',                s = cv_scorer.r2();
                case 'pearson_r',         s = cv_scorer.pearson_r();
                otherwise
                    error('cv_scorer:UnknownName', ...
                        'No scorer named ''%s''. Available: %s', ...
                        name, strjoin(cv_scorer.names(), ', '));
            end
        end

        function L = names()
            L = {'accuracy','balanced_accuracy','f1','roc_auc','log_loss', ...
                 'mse','rmse','mae','r2','pearson_r'};
        end

        % -------------------- Default-by-task --------------------
        function s = default_for_task(task)
            switch lower(task)
                case 'classification', s = cv_scorer.balanced_accuracy();
                case 'regression',     s = cv_scorer.r2();
                otherwise
                    error('cv_scorer:UnknownTask', 'Unknown task: %s', task);
            end
        end

        % -------------------- Metric implementations --------------------
        function v = fn_balanced_accuracy(Y, yh, ~)
            classes = unique(Y);
            recalls = zeros(numel(classes), 1);
            for i = 1:numel(classes)
                msk = Y == classes(i);
                if any(msk)
                    recalls(i) = mean(yh(msk) == classes(i));
                end
            end
            v = mean(recalls);
        end

        function v = fn_f1_binary(Y, yh, ~)
            classes = unique(Y);
            if numel(classes) ~= 2
                error('cv_scorer:f1:Binary', 'f1 supports binary classification only.');
            end
            pos = max(classes); % convention: larger class is positive
            tp = sum(yh == pos & Y == pos);
            fp = sum(yh == pos & Y ~= pos);
            fn = sum(yh ~= pos & Y == pos);
            if tp + fp == 0 || tp + fn == 0
                v = 0; return
            end
            prec = tp / (tp + fp);
            rec  = tp / (tp + fn);
            if prec + rec == 0
                v = 0;
            else
                v = 2 * prec * rec / (prec + rec);
            end
        end

        function v = fn_roc_auc(Y, ~, scores)
            if isempty(scores)
                error('cv_scorer:roc_auc:NeedsScores', ...
                      'roc_auc requires the third (continuous scores) argument.');
            end
            scores = scores(:); Y = Y(:);
            classes = unique(Y);
            if numel(classes) ~= 2
                error('cv_scorer:roc_auc:Binary', 'roc_auc supports binary classification only.');
            end
            pos = max(classes);
            pos_msk = Y == pos;
            if ~any(pos_msk) || all(pos_msk)
                v = 0.5; return
            end
            % AUC via Mann-Whitney U / Wilcoxon
            ranks = tiedrank(scores);
            n_pos = sum(pos_msk);
            n_neg = numel(Y) - n_pos;
            sum_ranks_pos = sum(ranks(pos_msk));
            v = (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg);
        end

        function v = fn_log_loss(Y, ~, scores)
            if isempty(scores)
                error('cv_scorer:log_loss:NeedsScores', ...
                      'log_loss requires the third (continuous probability) argument.');
            end
            scores = scores(:); Y = Y(:);
            p = max(min(scores, 1 - 1e-15), 1e-15);  % clip
            classes = unique(Y); pos = max(classes);
            y_bin = (Y == pos);
            v = -mean(y_bin .* log(p) + (1 - y_bin) .* log(1 - p));
        end

        function v = fn_r2(Y, yh, ~)
            ss_res = sum((Y - yh).^2);
            ss_tot = sum((Y - mean(Y)).^2);
            if ss_tot == 0
                v = NaN;
            else
                v = 1 - ss_res / ss_tot;
            end
        end

        function r = safe_corr(Y, yh)
            if numel(Y) < 2 || std(Y) == 0 || std(yh) == 0
                r = NaN; return
            end
            r = corr(Y(:), yh(:), 'Rows', 'complete');
        end
    end

end
