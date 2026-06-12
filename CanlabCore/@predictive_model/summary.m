function s = summary(obj, varargin)
% summary  Human-readable overview of a predictive_model: how it was built,
% what inference is available, and how well it performs.
%
% Prints (and optionally returns) a compact report covering:
%   - identity:    algorithm, task, fit provenance (insample / crossval / test)
%   - data:        number of observations, features, and grouping/within-person
%   - validation:  cross-validation scheme (splitter type, folds, scorer)
%   - inference:   whether bootstrap weights, a permutation test, calibration,
%                  or cross-classification results are present
%   - performance: the model-type-relevant metrics block from report_accuracy
%                  (classification: accuracy / balanced acc / AUC / sens /
%                  spec / PPV / NPV / d; regression: r / predicted R^2 /
%                  out-of-sample R^2 / RMSE / MAE)
%
% :Usage:
% ::
%     summary(pm);              % print the report
%     s = summary(pm);          % struct: s.provenance + s.accuracy (still prints)
%     s = summary(pm, 'noverbose');   % struct only, no printing
%
% :Inputs:
%
%   **obj:**
%        a @predictive_model in any state (more is reported as more is filled
%        in: fit -> crossval -> bootstrap / permutation_test / calibrate).
%
% :Optional Inputs (name/value):
%
%   **'doverbose' / 'verbose' / 'noverbose':**
%        control printing (default print on). 'noverbose' suppresses output.
%
% :Outputs:
%
%   **s:**
%        struct with fields .provenance (algorithm, task, fit_type, n_obs,
%        n_features, n_groups, cv_type, nfolds, scorer, has_bootstrap,
%        has_permutation, has_calibration, has_cross_classify) and .accuracy
%        (the struct from report_accuracy).
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm', 'noverbose');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm  = predictive_model('algorithm','svm','task','classification');
%     pm  = crossval(pm, X, Y, 'groups', id);
%     pm  = bootstrap(pm, X, Y, 'nboot', 1000, 'groups', id);
%     summary(pm);
%
% :See also:
%   report_accuracy, crossval, bootstrap, permutation_test, predictive_model

    doverbose = parse_verbose(varargin{:});

    prov = local_provenance(obj);

    if doverbose
        local_print(obj, prov);
    end

    % Compute the metrics once; report_accuracy prints the performance block
    % when doverbose (kept after the provenance header).
    acc = report_accuracy(obj, 'doverbose', doverbose);

    if nargout > 0
        s = struct('provenance', prov, 'accuracy', acc);
    end
end


% -------------------------------------------------------------------------
function prov = local_provenance(obj)
    prov.algorithm = char(string(obj.algorithm));
    prov.task      = char(string(obj.task));
    prov.fit_type  = char(string(obj.fit_type));

    prov.n_obs = 0;
    if ~isempty(obj.Y), prov.n_obs = sum(~isnan(obj.Y(:))); end

    prov.n_features = NaN;
    if isstruct(obj.weights) && isfield(obj.weights, 'w') && ~isempty(obj.weights.w)
        prov.n_features = numel(obj.weights.w);
    elseif ~isempty(obj.omitted_features)
        prov.n_features = numel(obj.omitted_features);
    end

    prov.n_groups = NaN;
    if ~isempty(obj.id)
        prov.n_groups = numel(unique(obj.id(~ismissing_safe(obj.id))));
    end
    prov.within_person = isstruct(obj.diagnostics) ...
        && isfield(obj.diagnostics, 'mult_obs_within_person') ...
        && isequal(obj.diagnostics.mult_obs_within_person, true);

    % Cross-validation scheme.
    prov.cv_type = '';
    if ~isempty(obj.cv) && isprop_or_field(obj.cv, 'type')
        prov.cv_type = char(string(obj.cv.type));
    end
    prov.nfolds = NaN;
    if isstruct(obj.cv_partition) && isfield(obj.cv_partition, 'nfolds') ...
            && ~isempty(obj.cv_partition.nfolds)
        prov.nfolds = obj.cv_partition.nfolds;
    end
    prov.scorer = '';
    if ~isempty(obj.scorer) && isprop_or_field(obj.scorer, 'name')
        prov.scorer = char(string(obj.scorer.name));
    end

    % Inference availability.
    prov.has_bootstrap = isstruct(obj.weights) && isfield(obj.weights, 'p') ...
        && ~isempty(obj.weights.p);
    prov.has_permutation = isstruct(obj.permutation_results) ...
        && ~isempty(fieldnames(obj.permutation_results));
    prov.has_calibration = isstruct(obj.diagnostics) ...
        && isfield(obj.diagnostics, 'calibration') ...
        && ~isempty(obj.diagnostics.calibration);
    prov.has_cross_classify = isstruct(obj.cross_classify) ...
        && ~isempty(fieldnames(obj.cross_classify));
end


% -------------------------------------------------------------------------
function local_print(obj, prov) %#ok<INUSL>
    fprintf('\n=== predictive_model summary ===\n');
    fprintf('  Algorithm : %s   (task: %s)\n', prov.algorithm, prov.task);

    switch lower(prov.fit_type)
        case 'crossval', ft = 'cross-validated';
        case 'insample', ft = 'in-sample (full fit, no CV)';
        case 'test',     ft = 'applied to held-out test data';
        otherwise,       ft = prov.fit_type;
    end
    if isempty(ft), ft = '(not yet fit)'; end
    fprintf('  Fit       : %s\n', ft);

    feat = '?'; if ~isnan(prov.n_features), feat = sprintf('%d', prov.n_features); end
    fprintf('  Data      : %d observations, %s features', prov.n_obs, feat);
    if ~isnan(prov.n_groups)
        fprintf(', %d groups%s', prov.n_groups, ...
            ternary(prov.within_person, ' (within-person design)', ''));
    end
    fprintf('\n');

    if ~isempty(prov.cv_type) || ~isnan(prov.nfolds)
        folds = '?'; if ~isnan(prov.nfolds), folds = sprintf('%d', prov.nfolds); end
        sc = ''; if ~isempty(prov.scorer), sc = sprintf(', scorer=%s', prov.scorer); end
        fprintf('  CV scheme : %s, %s folds%s\n', ...
            ternary(isempty(prov.cv_type), 'kfold', prov.cv_type), folds, sc);
    end

    avail = {};
    if prov.has_bootstrap,      avail{end+1} = 'bootstrap weights'; end
    if prov.has_permutation,    avail{end+1} = 'permutation test'; end
    if prov.has_calibration,    avail{end+1} = 'calibration'; end
    if prov.has_cross_classify, avail{end+1} = 'cross-classification'; end
    if isempty(avail)
        fprintf('  Inference : none yet (run bootstrap / permutation_test / calibrate)\n');
    else
        fprintf('  Inference : %s\n', strjoin(avail, ', '));
    end
end


% -------------------------------------------------------------------------
function tf = isprop_or_field(x, name)
    tf = (isobject(x) && isprop(x, name)) || (isstruct(x) && isfield(x, name));
end

function tf = ismissing_safe(x)
    try
        tf = ismissing(x);
    catch
        tf = false(size(x));
    end
end

function out = ternary(tf, a, b)
    if tf, out = a; else, out = b; end
end

function doverbose = parse_verbose(varargin)
    doverbose = true;
    for i = 1:numel(varargin)
        a = varargin{i};
        if (ischar(a) || isstring(a))
            if strcmpi(a, 'noverbose'), doverbose = false;
            elseif any(strcmpi(a, {'doverbose','verbose'})) && i < numel(varargin) ...
                    && (islogical(varargin{i+1}) || isnumeric(varargin{i+1}))
                doverbose = logical(varargin{i+1});
            end
        end
    end
end
