function [allcverr, allyhat, pm_all, results] = predict_test_suite(dat, varargin)
% Cross-validate a set of @predictive_model algorithms on an fmri_data
% object, plot predicted-vs-observed for each, and tabulate performance.
%
% :Usage:
% ::
%
%     [allcverr, allyhat, pm_all, results] = predict_test_suite(dat, [optional inputs])
%
% This is the new-API port of the legacy predict_test_suite. Instead of
% calling fmri_data.predict repeatedly, it constructs one
% @predictive_model per algorithm, cross-validates each over a shared
% cv_splitter (so the folds are identical across algorithms), and
% compares cross-validated score, error, and (optionally) the number of
% bootstrap-FDR-significant features.
%
% :Inputs:
%
%   **dat:**
%        an fmri_data object. dat.Y must be assigned with either binary
%        (classification) or continuous (regression) outcomes. Algorithms
%        are chosen automatically from the number of unique values in Y.
%
% :Optional Inputs:
%
%   **'nfolds':**
%        Followed by either a scalar k (k-fold cross-validation, default 5)
%        or a vector of integer fold/holdout ids (one per observation, e.g.
%        subject id) used verbatim as a custom partition.
%
%   **'algorithms':**
%        Followed by a cell array of @predictive_model registry names to
%        override the automatic selection (e.g. {'svm','logistic'}).
%
%   **'bootstrap':**
%        Followed by a number of bootstrap samples. When given, each model
%        is bootstrapped and the table reports the count of FDR-significant
%        features (weights.fdr_sig). Off by default (bootstrapping is slow).
%
%   **'quick':**
%        Skip the extended hyperparameter-sweep panel.
%
%   **'noplot'** / **'doplot':**
%        Suppress / force plotting (default: plot).
%
% :Outputs:
%
%   **allcverr:**
%        [1 x n_methods] cross-validated error (misclassification rate for
%        classification, mean squared error for regression).
%
%   **allyhat:**
%        [n_obs x n_methods] cross-validated predictions per algorithm.
%
%   **pm_all:**
%        {1 x n_methods} cell of the cross-validated @predictive_model
%        objects (and bootstrapped, if 'bootstrap' was requested).
%
%   **results:**
%        a table with one row per algorithm: cv_score, cv_error, and
%        n_fdr_sig (NaN unless 'bootstrap' was requested).
%
% :Examples:
% ::
%
%    % Classification on a binary outcome
%    dat = load_image_set('DPSP_hotwarm');           % fmri_data, binary Y
%    [cverr, yhat, pm_all, results] = predict_test_suite(dat, 'nfolds', dat.metadata_table.subject_id);
%    disp(results)
%
%    % Regression, with bootstrap-FDR feature counts
%    [cverr, yhat, pm_all, results] = predict_test_suite(dat, 'bootstrap', 1000);
%
% :See also:
%   predictive_model, crossval, bootstrap, fmri_data.predict
%
% ..
%    Tor Wager. Original version Nov 2012; ported to the @predictive_model
%    API in the predictive-model object-development branch.
% ..

% ---------------------------------------------------------------------
% Parse inputs
% ---------------------------------------------------------------------
p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'nfolds',     5);
addParameter(p, 'algorithms', {});
addParameter(p, 'bootstrap',  0);
addParameter(p, 'quick',      false);
addParameter(p, 'doplot',     true);

% Accept bare flags ('quick', 'noplot', 'doplot') in CANlab style.
keep = {};
for i = 1:numel(varargin)
    if ischar(varargin{i}) || isstring(varargin{i})
        switch lower(char(varargin{i}))
            case 'quick',  varargin{i} = []; keep = [keep, {'quick', true}];   %#ok<AGROW>
            case 'noplot', varargin{i} = []; keep = [keep, {'doplot', false}]; %#ok<AGROW>
            case 'doplot', varargin{i} = []; keep = [keep, {'doplot', true}];  %#ok<AGROW>
        end
    end
end
varargin = varargin(~cellfun(@isempty, varargin));
parse(p, varargin{:}, keep{:});

nfolds     = p.Results.nfolds;
algorithms = p.Results.algorithms;
nboot      = p.Results.bootstrap;
do_quick   = logical(p.Results.quick);
do_plot    = logical(p.Results.doplot);

% ---------------------------------------------------------------------
% Extract X, Y; pick task and algorithms
% ---------------------------------------------------------------------
X = double(dat.dat');
Y = double(dat.Y(:));

nlevels = numel(unique(Y(~isnan(Y))));

switch nlevels
    case {0, 1}
        error('predict_test_suite:NoVariance', 'No variance in dat.Y.');
    case 2
        task = 'classification';
        if isempty(algorithms), algorithms = {'svm', 'linear_svm', 'logistic'}; end
        scoring = 'accuracy';
    otherwise
        task = 'regression';
        if isempty(algorithms), algorithms = {'svr', 'linear_svr', 'lasso', 'ridge'}; end
        scoring = '';   % default-for-task (r2)
end

% ---------------------------------------------------------------------
% Shared cross-validation splitter (identical folds for every algorithm)
% ---------------------------------------------------------------------
if isscalar(nfolds)
    if strcmp(task, 'classification')
        cv = cv_splitter.stratified_kfold(nfolds);
    else
        cv = cv_splitter.kfold(nfolds);
    end
    cv_descrip = sprintf('%d-fold', nfolds);
else
    cv = cv_splitter.custom_partition(nfolds(:));
    cv_descrip = sprintf('%d custom folds', numel(unique(nfolds)));
end

% ---------------------------------------------------------------------
% Cross-validate each algorithm
% ---------------------------------------------------------------------
n_methods = numel(algorithms);
pm_all    = cell(1, n_methods);
allyhat   = nan(numel(Y), n_methods);
allcverr  = nan(1, n_methods);
cv_score  = nan(1, n_methods);
n_fdr_sig = nan(1, n_methods);

if do_plot
    create_figure('predicted-actual corr', 2, max(4, n_methods + 1));
end

for i = 1:n_methods

    try
        pm = predictive_model('algorithm', algorithms{i}, 'task', task);

        if isempty(scoring)
            pm = crossval(pm, X, Y, 'cv', cv);
        else
            pm = crossval(pm, X, Y, 'cv', cv, 'scoring', scoring);
        end

        yhat = pm.fitted_values.yfit;
        allyhat(:, i) = yhat;
        valid = ~isnan(yhat) & ~isnan(Y);

        % cv_score = the scorer's mean across folds.
        cv_score(i) = predictive_model.metric_value(pm.error_metrics, pm.scorer.name);

        % cv_error: mcr (classification) or mse (regression).
        if strcmp(task, 'classification')
            allcverr(i) = mean(Y(valid) ~= yhat(valid));
        else
            allcverr(i) = mean((Y(valid) - yhat(valid)).^2);
        end

        if nboot > 0
            pm = bootstrap(pm, X, Y, 'nboot', nboot);
            if isfield(pm.weights, 'fdr_sig') && ~isempty(pm.weights.fdr_sig)
                n_fdr_sig(i) = sum(pm.weights.fdr_sig);
            end
        end

        pm_all{i} = pm;

        if do_plot
            subplot(2, max(4, n_methods + 1), i);
            local_scatter(yhat, Y, task);
            title(algorithms{i}, 'Interpreter', 'none');
            local_plot_folds(yhat, Y, pm.cv_partition);
            drawnow
        end

    catch ME
        fprintf('ERROR running %s: %s\n', algorithms{i}, ME.message);
    end

end

% ---------------------------------------------------------------------
% Summary table + cross-method prediction correlations
% ---------------------------------------------------------------------
results = table(algorithms(:), cv_score(:), allcverr(:), n_fdr_sig(:), ...
    'VariableNames', {'algorithm', 'cv_score', 'cv_error', 'n_fdr_sig'});

fprintf('\nCross-validated prediction (%s, %s):\n', task, cv_descrip);
disp(results);

ok = all(~isnan(allyhat), 1);
if sum(ok) > 1
    rmtx = corrcoef(allyhat(:, ok));
    disp('Correlations among predicted values across methods:');
    print_matrix(rmtx, algorithms(ok), algorithms(ok));
end

if do_quick || ~do_plot
    return
end

% ---------------------------------------------------------------------
% Extended panel: regularization (Lambda) sweep on a representative
% penalized algorithm.
% ---------------------------------------------------------------------
try
    if strcmp(task, 'classification')
        sweep_alg = 'logistic';
    else
        sweep_alg = 'lasso';
    end

    lambdas = logspace(-4, 0, 8);
    sweep_score = nan(size(lambdas));
    for k = 1:numel(lambdas)
        pmk = predictive_model('algorithm', sweep_alg, 'task', task, ...
            'modeloptions', {'Lambda', lambdas(k)});
        if isempty(scoring)
            pmk = crossval(pmk, X, Y, 'cv', cv);
        else
            pmk = crossval(pmk, X, Y, 'cv', cv, 'scoring', scoring);
        end
        sweep_score(k) = predictive_model.metric_value(pmk.error_metrics, pmk.scorer.name);
    end

    subplot(2, max(4, n_methods + 1), n_methods + 1);
    semilogx(lambdas, sweep_score, 'o-', 'Color', [.3 .3 .3], ...
        'MarkerFaceColor', [.5 .5 1], 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Lambda (regularization)');
    if strcmp(task, 'classification')
        ylabel('CV accuracy');
    else
        ylabel('CV R^2');
    end
    title(sprintf('Reg. sweep: %s', sweep_alg), 'Interpreter', 'none');

catch ME
    fprintf('Regularization sweep skipped: %s\n', ME.message);
end

end % main function


% =====================================================================
% local helpers
% =====================================================================
function local_scatter(yhat, Y, task)
    valid = ~isnan(yhat) & ~isnan(Y);
    if strcmp(task, 'classification')
        classes = unique(Y(valid));
        colors  = scn_standard_colors(numel(classes));
        hold on
        for c = 1:numel(classes)
            msk = valid & Y == classes(c);
            plot(yhat(msk) + 0.05*randn(sum(msk),1), Y(msk), 'o', ...
                'MarkerFaceColor', colors{c}, 'MarkerEdgeColor', 'none');
        end
        xlabel('Predicted class (cv)'); ylabel('Observed class');
    else
        plot_correlation_samefig(yhat(valid), Y(valid));
        xlabel('Predicted (cv)'); ylabel('Observed');
    end
end


function local_plot_folds(yhat, Y, cv_partition)
% Color-code held-out predictions by fold (regression only).
    if ~isfield(cv_partition, 'teIdx') || isempty(cv_partition.teIdx), return; end
    if numel(unique(Y(~isnan(Y)))) <= 2, return; end   % skip for classification

    teIdx  = cv_partition.teIdx;
    colors = scn_standard_colors(numel(teIdx));
    hold on
    for i = 1:numel(teIdx)
        wh = find(teIdx{i});
        if numel(wh) < 3, continue; end
        yf = yhat(wh); y = Y(wh);
        ok = ~isnan(yf) & ~isnan(y);
        if sum(ok) < 3, continue; end
        plot(yf(ok), y(ok), 'o', 'MarkerFaceColor', colors{i});
        b = glmfit(yf(ok), y(ok));
        plot([min(yf(ok)) max(yf(ok))], ...
             [b(1)+b(2)*min(yf(ok)) b(1)+b(2)*max(yf(ok))], ...
             'Color', colors{i}, 'LineWidth', 2);
    end
    axis tight
end
