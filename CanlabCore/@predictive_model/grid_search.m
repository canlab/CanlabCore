function obj = grid_search(obj, X, Y, param_grid, varargin)
% grid_search  Exhaustive hyperparameter search via cross-validation.
%
% :Usage:
% ::
%     param_grid.BoxConstraint    = [0.1 1 10 100];
%     param_grid.KernelScale      = [0.5 1 2 5];
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = grid_search(pm, X, Y, param_grid);
%
%     % nested CV: grid_search inside crossval
%     pm = predictive_model('algorithm','svr','task','regression');
%     pm = crossval(pm, X, Y, 'cv', cv_splitter.kfold(5));
%     % ... after which pm.modeloptions could be updated by grid_search
%     %     on each outer fold (not yet automated — see plan §G6).
%
% Pipeline:
%   1. Build the Cartesian product of param_grid values.
%   2. For each combination:
%        - clone(obj); append the (name, value) pair to modeloptions
%        - crossval on (X, Y) with obj.cv (or default)
%        - record mean CV score
%   3. Pick the combination with the best score (argmax if
%      scorer.greater_is_better, else argmin).
%   4. Refit obj with the best modeloptions via crossval, so
%      obj.error_metrics holds the chosen model's CV performance.
%   5. Stash the full grid scores in obj.diagnostics.grid_search.
%
% :Inputs:
%   param_grid   struct mapping option-name -> vector (or cell) of values.
%                Each value is forwarded to the MATLAB fitter as a
%                modeloptions name/value pair.
%
% :Optional Inputs (name/value):
%   'groups'     grouping vector for grouped CV
%   'verbose'    default true
%
% :Outputs:
%
%   **obj:**
%        the @predictive_model refit with the best hyperparameters:
%        obj.modeloptions extended with the winning name/value pairs,
%        obj.error_metrics holding the chosen model's CV performance, and
%        obj.diagnostics.grid_search a struct with .param_names, .combos,
%        .scores, .best_idx, .best_score.
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm = predictive_model('algorithm','svm','task','classification');
%     grid.BoxConstraint = [0.1 1 10];
%     pm = grid_search(pm, X, Y, grid, 'groups', id);
%     pm.diagnostics.grid_search.best_score
%
% :See also:
%   crossval, fit, stability_selection, cv_splitter

    pi = inputParser; pi.KeepUnmatched = true;
    addParameter(pi, 'groups',  []);
    addParameter(pi, 'verbose', true);
    parse(pi, varargin{:});

    groups  = pi.Results.groups;
    verbose = pi.Results.verbose;

    names = fieldnames(param_grid);
    if isempty(names)
        error('predictive_model:grid_search:EmptyGrid', ...
            'param_grid must contain at least one parameter.');
    end

    % Build vectors of values per param.
    vals_per_name = cell(numel(names), 1);
    for i = 1:numel(names)
        v = param_grid.(names{i});
        if isnumeric(v) || islogical(v)
            vals_per_name{i} = num2cell(v(:));
        elseif iscell(v)
            vals_per_name{i} = v(:);
        else
            error('predictive_model:grid_search:BadValue', ...
                'param_grid.%s must be numeric, logical, or cell.', names{i});
        end
    end

    % Cartesian product via ndgrid.
    grid_cells = cell(numel(names), 1);
    [grid_cells{:}] = ndgrid(vals_per_name{:});
    n_combos = numel(grid_cells{1});

    % Pick scorer to determine argmax/min direction.
    if isempty(obj.scorer)
        if isempty(obj.task)
            % infer from Y
            if numel(unique(Y(~isnan(Y)))) <= 2
                obj.task = 'classification';
            else
                obj.task = 'regression';
            end
        end
        obj.scorer = cv_scorer.default_for_task(obj.task);
    end

    scores = nan(n_combos, 1);

    if verbose
        fprintf('grid_search: %d combinations of %d params...\n', n_combos, numel(names));
    end

    for k = 1:n_combos
        m = clone(obj);
        extra_opts = cell(1, 2 * numel(names));
        for j = 1:numel(names)
            cell_jk = grid_cells{j};  % ndgrid cell array of param j
            v = cell_jk{k};
            if iscell(v), v = v{1}; end
            extra_opts{2*j-1} = names{j};
            extra_opts{2*j}   = v;
        end
        m.modeloptions = [m.modeloptions, extra_opts];

        try
            m = crossval(m, X, Y, 'groups', groups);
            scores(k) = m.error_metrics.(m.scorer.name).value;
            if verbose
                fprintf('  [%d/%d] %s = %.4f\n', k, n_combos, ...
                    grid_combo_string(names, extra_opts), scores(k));
            end
        catch ME
            if verbose
                fprintf('  [%d/%d] FAILED: %s\n', k, n_combos, ME.message);
            end
        end
    end

    % Pick best.
    if obj.scorer.greater_is_better
        [best_score, best_idx] = max(scores);
    else
        [best_score, best_idx] = min(scores);
    end

    if isnan(best_score)
        error('predictive_model:grid_search:AllFailed', ...
            'All %d grid combinations failed.', n_combos);
    end

    % Apply best params and re-run crossval for final populated obj.
    best_opts = cell(1, 2 * numel(names));
    for j = 1:numel(names)
        cell_jk = grid_cells{j};
        v = cell_jk{best_idx};
        if iscell(v), v = v{1}; end
        best_opts{2*j-1} = names{j};
        best_opts{2*j}   = v;
    end
    obj.modeloptions = [obj.modeloptions, best_opts];
    obj = crossval(obj, X, Y, 'groups', groups);

    % Stash grid results in diagnostics.
    obj.diagnostics.grid_search = struct( ...
        'param_names', {names}, ...
        'scores',      scores, ...
        'best_idx',    best_idx, ...
        'best_score',  best_score, ...
        'best_params', {best_opts});

    obj.history{end+1, 1} = sprintf('grid_search: %d combos, best %s = %.4f', ...
        n_combos, obj.scorer.name, best_score);

    if verbose
        fprintf('grid_search: best %s = %.4f at %s\n', ...
            obj.scorer.name, best_score, grid_combo_string(names, best_opts));
    end
end


% --- helper ---------------------------------------------------------------
function s = grid_combo_string(names, opts)
    parts = cell(1, numel(names));
    for j = 1:numel(names)
        v = opts{2*j};
        if isnumeric(v)
            parts{j} = sprintf('%s=%g', names{j}, v);
        elseif ischar(v) || isstring(v)
            parts{j} = sprintf('%s=%s', names{j}, char(v));
        else
            parts{j} = sprintf('%s=<%s>', names{j}, class(v));
        end
    end
    s = strjoin(parts, ', ');
end
