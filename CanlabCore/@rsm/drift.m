function varargout = drift(obj, varargin)
% drift  Quantify within-condition pattern stability across replicates.
%
% Generalizes the `08192024 Representational Drift in RSA` workflow: for each
% condition (or grouping), pulls the per-replicate self-similarity / row from
% the rsm and quantifies how it changes across replicates (sessions/runs).
%
% Two operating modes
% -------------------
%   'reference','first' (default) -- correlate each replicate's RSM row with
%                                    the first replicate's. Returns a per-
%                                    replicate-per-condition matrix of
%                                    correlations.
%
%   'reference','mean'           -- correlate each replicate's RSM with the
%                                   mean RSM across the other replicates
%                                   (leave-one-out).
%
% Optional linear fit
% -------------------
%   'fit','linear'  -- additionally fit a per-condition linear model
%                      similarity ~ replicate_index. Returns slopes, intercepts,
%                      R-squared, p-values in a side table.
%
% Usage
% -----
%   % Per-condition similarity-to-first-replicate
%   D = R.drift();                     % returns table: replicate | condition | similarity
%
%   % With replicate_table info (subject, session)
%   D = R.drift('reference','first')
%
%   % Linear fit of similarity over replicate index, per condition
%   [D, fit] = R.drift('fit','linear');
%
% Inputs
% ------
%   varargin:
%     'reference'  'first' (default) | 'mean'
%     'fit'        'none' (default) | 'linear'
%     'metric'     'pearson' (default) | 'spearman'
%
% Output
% ------
%   If 'fit' = 'none':
%       table with columns: replicate_index | condition | similarity (+ replicate_table cols)
%   If 'fit' = 'linear':
%       [drift_table, fit_table] where fit_table is:
%         condition | slope | intercept | R_squared | P_value

if builtin('numel', obj) > 1
    error('rsm:drift:nonScalar', 'expects a scalar rsm.');
end

p = inputParser;
p.addParameter('reference', 'first',   @(x) ischar(x) || isstring(x));
p.addParameter('fit',       'none',    @(x) ischar(x) || isstring(x));
p.addParameter('metric',    'pearson', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

[k, ~, N] = size(obj.dat);
if N < 2
    error('rsm:drift:tooFewReplicates', 'Need at least 2 replicates; got %d.', N);
end

% For each condition (row), build the "drift" series: similarity of that
% condition's relationships to the reference, across replicates.
% Per condition c: at each replicate n, take RSM(c, others, n) and correlate
% with the reference RSM(c, others, *).

ref_mode = lower(char(opt.reference));
metric = lower(char(opt.metric));
corr_type = 'Pearson'; if strcmp(metric, 'spearman'), corr_type = 'Spearman'; end

cond_labels = obj.labels;
if isempty(cond_labels)
    cond_labels = arrayfun(@(i) sprintf('cond_%d', i), (1:k)', 'UniformOutput', false);
end

% Pre-compute reference row vectors per condition
ref_rows = nan(k, k - 1);   % row per condition, excluding self
for c = 1:k
    other = [1:c-1, c+1:k];
    switch ref_mode
        case 'first'
            ref_rows(c, :) = obj.dat(c, other, 1);
        case 'mean'
            ref_rows(c, :) = mean(obj.dat(c, other, :), 3, 'omitnan');
        otherwise
            error('rsm:drift:badReference', ...
                'reference must be ''first'' or ''mean''.');
    end
end

% Build the per-(condition, replicate) similarity series
n_rows = N * k;
rep_idx = zeros(n_rows, 1);
cond_col = cell(n_rows, 1);
sim_col = nan(n_rows, 1);
row = 0;
for n = 1:N
    for c = 1:k
        other = [1:c-1, c+1:k];
        v = squeeze(obj.dat(c, other, n))';
        row = row + 1;
        rep_idx(row) = n;
        cond_col{row} = cond_labels{c};
        if strcmp(ref_mode, 'mean')
            % Leave-one-out: exclude self-replicate from the reference
            ref = mean(obj.dat(c, other, setdiff(1:N, n)), 3, 'omitnan');
            sim_col(row) = nan_safe_corr(v, ref, corr_type);
        else
            sim_col(row) = nan_safe_corr(v, ref_rows(c, :)', corr_type);
        end
    end
end

drift_table = table(rep_idx, cond_col, sim_col, ...
    'VariableNames', {'replicate_index', 'condition', 'similarity'});

% Attach replicate_table columns
if ~isempty(obj.replicate_table)
    rep_cols = obj.replicate_table.Properties.VariableNames;
    for i = 1:numel(rep_cols)
        v = obj.replicate_table.(rep_cols{i});
        if isnumeric(v)
            drift_table.(rep_cols{i}) = v(rep_idx);
        elseif iscell(v)
            drift_table.(rep_cols{i}) = v(rep_idx);
        elseif iscategorical(v) || isstring(v)
            drift_table.(rep_cols{i}) = v(rep_idx);
        end
    end
end

if strcmpi(opt.fit, 'none')
    varargout{1} = drift_table;
    if nargout >= 2, varargout{2} = table.empty; end
    return
end

if ~strcmpi(opt.fit, 'linear')
    error('rsm:drift:badFit', 'fit must be ''none'' or ''linear''.');
end

% Per-condition linear fit of similarity ~ replicate_index
fit_rows = cell(k, 1);
for c = 1:k
    is_c = strcmp(drift_table.condition, cond_labels{c});
    x = drift_table.replicate_index(is_c);
    y = drift_table.similarity(is_c);
    mask = isfinite(y);
    if nnz(mask) < 3
        fit_rows{c} = table({cond_labels{c}}, NaN, NaN, NaN, NaN, ...
            'VariableNames', {'condition', 'slope', 'intercept', 'R_squared', 'P_value'});
        continue
    end
    mdl = fitlm(x(mask), y(mask));
    fit_rows{c} = table({cond_labels{c}}, mdl.Coefficients.Estimate(2), ...
                       mdl.Coefficients.Estimate(1), mdl.Rsquared.Ordinary, ...
                       mdl.Coefficients.pValue(2), ...
        'VariableNames', {'condition', 'slope', 'intercept', 'R_squared', 'P_value'});
end
fit_table = vertcat(fit_rows{:});

varargout{1} = drift_table;
if nargout >= 2
    varargout{2} = fit_table;
end

end


function r = nan_safe_corr(a, b, corr_type)
% Pairwise-complete correlation with a minimum-overlap floor of 3 finite
% pairs. Avoids NaN propagation when one input has missing entries (e.g. a
% session with one missing condition produces a partial NaN row).
a = a(:); b = b(:);
m = isfinite(a) & isfinite(b);
if nnz(m) < 3
    r = NaN; return
end
r = corr(a(m), b(m), 'Type', corr_type);
end
