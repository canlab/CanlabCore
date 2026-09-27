function T = reliability_by_grouping(obj, varargin)
% reliability_by_grouping  ICC reliability for each grouping (e.g. per superordinate condition).
%
% For each grouping defined in R.groupings (or specified explicitly), restricts
% the RSM to that subset of rows/cols and computes ICC across the replicate axis.
% Reproduces the per-condition / per-bodysite loops in `01282025 RSM Reliability`.
%
% Usage
% -----
%   T = R.reliability_by_grouping()                                % uses R.groupings
%   T = R.reliability_by_grouping('groupings',{'hot','warm','imagine'})
%   T = R.reliability_by_grouping('icc_type','2-k', 'by','subject_id')
%
% Optional name-value
% -------------------
%   'groupings'  cellstr of grouping names (default: all in R.groupings)
%   'icc_type'   '3-k' (default) | '2-k' | etc. (see reliability())
%   'by'         '' (default) or replicate_table column for per-(group, by) ICC
%   'transform'  'auto' (default) | 'fisherz' | 'none'
%
% Output
% ------
%   T  table:
%        - if 'by' is empty: {Grouping, ICC, n_replicates, n_cells}
%        - if 'by' is set:   {Grouping, <by_col>, ICC, n_replicates, n_cells}

if builtin('numel', obj) > 1
    error('rsm:reliability_by_grouping:nonScalar', 'expects a scalar rsm.');
end

p = inputParser;
p.addParameter('groupings',          {},     @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('icc_type',           '3-k',  @(x) ischar(x) || isstring(x));
p.addParameter('by',                 '',     @(x) ischar(x) || isstring(x));
p.addParameter('transform',          'auto', @(x) ischar(x) || isstring(x));
p.addParameter('replicate_nan_frac', 0.5,    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('pool',               'auto', @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'auto','subject','replicate'}));
p.parse(varargin{:});
opt = p.Results;

% Resolve which groupings to use
gp = cellstr(opt.groupings);
if isempty(gp)
    if isempty(obj.groupings) || isempty(fieldnames(obj.groupings))
        error('rsm:reliability_by_grouping:noGroupings', ...
            'R.groupings is empty; pass ''groupings'' explicitly or attach R.groupings first.');
    end
    gp = fieldnames(obj.groupings);
end

n_g = numel(gp);
rows = cell(n_g, 1);

% Consolidate per-call warnings into one summary at the end. Save current
% warning state, suppress the dropped-replicate spam, then restore.
warn_ids = {'rsm:reliability:droppedReplicates', ...
            'rsm:reliability:tooFewReplicates', ...
            'rsm:reliability:tooManyNaN'};
saved_states = cell(numel(warn_ids), 1);
for i = 1:numel(warn_ids)
    saved_states{i} = warning('query', warn_ids{i});
    warning('off', warn_ids{i});
end
total_dropped_replicates = 0;
groupings_with_drops     = {};
groupings_with_small_n   = {};

for i = 1:n_g
    sub_R = obj.subset(gp{i});
    if size(sub_R, 1) < 2
        warning('rsm:reliability_by_grouping:singleton', ...
            'Grouping "%s" has fewer than 2 rows; skipping.', gp{i});
        continue
    end

    % Track small-n-cells groupings (consolidated warning at the end)
    n_cells_this = nchoosek(size(sub_R, 1), 2);
    if n_cells_this < 5
        groupings_with_small_n{end+1} = sprintf('%s(%d cells, k=%d)', ...
            gp{i}, n_cells_this, size(sub_R, 1)); %#ok<AGROW>
    end

    % Track dropped replicates by inspecting NaN content directly (cheap)
    [k_sub, ~, N_sub] = size(sub_R.dat);
    upper_mask = triu(true(k_sub), 1);
    nan_frac = zeros(N_sub, 1);
    for n = 1:N_sub
        slice = sub_R.dat(:,:,n);
        nan_frac(n) = mean(~isfinite(slice(upper_mask)));
    end
    n_drop = sum(nan_frac > opt.replicate_nan_frac);
    if n_drop > 0
        total_dropped_replicates = total_dropped_replicates + n_drop;
        groupings_with_drops{end+1} = sprintf('%s(%d/%d)', gp{i}, n_drop, N_sub); %#ok<AGROW>
    end

    inner = sub_R.reliability('icc_type', opt.icc_type, 'by', opt.by, ...
                              'transform', opt.transform, ...
                              'replicate_nan_frac', opt.replicate_nan_frac, ...
                              'pool', opt.pool);

    if isstruct(inner) && isfield(inner, 'per_subject')
        s = inner.summary;
        rows{i} = table({gp{i}}, s.mean, s.median, s.std, s.n_subjects, ...
            size(sub_R, 3), n_cells_this, ...
            'VariableNames', {'Grouping', 'Mean_ICC', 'Median_ICC', 'Std_ICC', ...
                              'n_subjects', 'n_replicates', 'n_cells'});
    elseif istable(inner)
        n_rows = height(inner);
        grp_col = repmat({gp{i}}, n_rows, 1);
        rows{i} = [table(grp_col, 'VariableNames', {'Grouping'}), inner];
    else
        rows{i} = table({gp{i}}, inner, size(sub_R,3), n_cells_this, ...
            'VariableNames', {'Grouping', 'ICC', 'n_replicates', 'n_cells'});
    end
end

% Restore prior warning states
for i = 1:numel(warn_ids)
    warning(saved_states{i}.state, warn_ids{i});
end

% Emit consolidated summaries
if ~isempty(groupings_with_small_n)
    warning('rsm:reliability_by_grouping:smallNcells', ...
        ['%d grouping(s) have <5 upper-triangle cells; ICC estimates are ', ...
         'statistically unreliable and can be unbounded. Affected: %s'], ...
        numel(groupings_with_small_n), ...
        strjoin(groupings_with_small_n, ', '));
end
if total_dropped_replicates > 0
    warning('rsm:reliability_by_grouping:droppedReplicates', ...
        ['Dropped replicates with > %.0f%% NaN cells across %d grouping(s): %s. ', ...
         'Total dropped (grouping x replicate) = %d. ', ...
         'Raise ''replicate_nan_frac'' to retain more, or use ', ...
         '''nan_policy'',''skip_replicate'' at compute_rsm time.'], ...
        100*opt.replicate_nan_frac, numel(groupings_with_drops), ...
        strjoin(groupings_with_drops, ', '), total_dropped_replicates);
end

T = vertcat(rows{:});
end
