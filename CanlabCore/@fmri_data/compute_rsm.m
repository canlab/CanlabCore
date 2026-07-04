function R = compute_rsm(dat, varargin)
% compute_rsm  Build an rsm object from an fmri_data object.
%
% Omnibus constructor that subsumes the four bespoke generate_RSA* wrappers
% in `C:\Temp\` and the inline RSM construction throughout the WASABI / Sun
% et al. workflow corpus. Drives all subsequent RSA inference machinery
% (cells, contrasts, ttest_contrasts, reliability, drift, rsa_lm, rsa_lme,
% compare_models, rsa_parcelwise, searchlight_rsa).
%
% Usage examples
% --------------
%   % One RSM per subject, conditions on the rows/cols, Spearman correlation
%   R = compute_rsm(dat, 'group_by','condition', 'subject_var','subject_id', ...
%                        'metric','spearman');
%
%   % Per-session 3D rsm (subjects x sessions stacked along dim 3)
%   R = compute_rsm(dat, 'level','session', 'group_by','condition', ...
%                        'subject_var','subject_id', 'session_var','session_number');
%
%   % Crossnobis with two-fold session split + session-difference whitening
%   R = compute_rsm(dat, 'metric','crossnobis', ...
%                        'group_by','condition', 'subject_var','subject_id', ...
%                        'fold_var','session_number', 'whiten','session_difference');
%
%   % Parcelwise: array of rsm objects, one per atlas parcel
%   R = compute_rsm(dat, 'parcellation', canlab2024_atlas, ...
%                        'group_by','condition', 'subject_var','subject_id');
%
%   % Sun et al. diagonal correction: replace same-bodysite diagonal with
%   % across-session same-bodysite mean
%   R = compute_rsm(dat, 'group_by','condition_x_bodysite', ...
%                        'subject_var','subject_id', ...
%                        'diagonal_correction','across_session_mean', ...
%                        'diagonal_group_by','bodysite', ...
%                        'session_var','session_number');
%
% Inputs
% ------
%   dat       fmri_data object. dat.metadata_table must contain all metadata
%             columns referenced by name in the optional arguments.
%
%   varargin name-value pairs:
%
%   Aggregation:
%     'group_by'             Metadata column whose unique values become the
%                            rows/cols of the RSM. Accepts either a single
%                            column name (string/char) or a cellstr of column
%                            names that get concatenated into a composite key.
%                            Example: {'condition','bodysite'} on Sun et al.
%                            data builds the canonical 24-row RSM (3 conditions
%                            x 8 bodysites).
%                            Default: '' (no grouping; one row per image -- image-level RSM).
%
%     'auto_groupings'       logical. When true (default), auto-builds R.groupings
%                            with one bare-name entry per unique value in each
%                            group_by column. For example, with
%                            group_by={'condition','bodysite'}, you get
%                            R.groupings.hot = 1:8, R.groupings.leftface = [1,9,17],
%                            etc. -- so reliability_by_grouping(), R.cells('hot','hot'),
%                            and R.contrast() work without manual setup.
%                            Set to false if you prefer to attach R.groupings yourself.
%     'condition_collapse'   'mean' | 'concat' | 'none'. How to aggregate when
%                            multiple images map to one condition. Default 'mean'.
%
%   Replicate axis:
%     'level'                'subject' (default) | 'session' | 'run' | 'collapsed' | 'image'
%                            How to split the data along the 3rd dim.
%                              'subject' -> one slice per subject
%                              'session' -> one slice per (subject, session)
%                              'run'     -> one slice per (subject, session, run)
%                              'collapsed' -> one slice (all data merged)
%                              'image'   -> no grouping; image-level RSM per subject
%     'subject_var'          Metadata column for subject IDs. Default 'subject_id'.
%     'session_var'          Metadata column for session number (required for
%                            level='session' or fold_var='session_number').
%                            Default 'session_number'.
%     'run_var'              Metadata column for run number (required for
%                            level='run').
%                            Default 'run_number'.
%
%   Metric:
%     'metric'               'correlation' (default) | 'spearman' | 'cosine' |
%                            'euclidean' | 'seuclidean' | 'mahalanobis' |
%                            'crossnobis' | 'cvcorr' | 'cvspearman'
%                            'cvcorr'/'cvspearman' are CROSS-VALIDATED
%                            correlation similarities: each cell (i,j) is the
%                            mean correlation between condition i and condition
%                            j taken from DIFFERENT folds (never the same fold).
%                            This removes within-fold shared structure (e.g.
%                            shared-run noise when several conditions co-occur
%                            in one run) while staying in correlation space --
%                            the diagonal is the cross-fold reliability, not 1.
%                            Requires 'fold_var'. Reproduces the Sun et al.
%                            (2026) BodyMap "exclude within-run correlations"
%                            RSM as a one-liner.
%
%   Cross-validated metrics (crossnobis / cvcorr / cvspearman):
%     'fold_var'             Metadata column defining cross-validation folds
%                            (required for these metrics). Default
%                            'session_number'. Special value 'occurrence'
%                            auto-builds folds by within-condition occurrence
%                            rank (each condition's 1st image -> fold 1, 2nd ->
%                            fold 2, ...), ordered by run_var when present.
%                            Use 'occurrence' when no metadata column lines the
%                            conditions up into shared folds (e.g. when run
%                            numbers differ across conditions).
%     'cv_scheme'            (cvcorr/cvspearman only) 'allpairs' (default) |
%                            'loo'. 'allpairs' correlates single-fold patterns
%                            for every ordered fold pair (most conservative,
%                            lowest power). 'loo' correlates each held-out fold
%                            against the mean of the other folds (higher SNR /
%                            MORE POWER, still never sharing a fold). Ignored by
%                            crossnobis (always all-pairs).
%
%   Whitening:
%     'whiten'               'none' (default) | 'within_subject' |
%                            'across_subject' | 'session_difference'
%                            within_subject and across_subject use covdiag.
%                            session_difference is the crossnobis whitening.
%     'whiten_method'        'covdiag' (default) | 'diag' | 'none'
%
%   Diagonal correction:
%     'diagonal_correction'  'none' (default) | 'across_session_mean' | 'nan'
%     'diagonal_group_by'    Metadata column to group by for diagonal correction.
%                            Default '' (uses 'group_by' if set).
%
%   Spatial restriction & preprocessing:
%     'mask'                 image_vector / atlas / region / char keyword.
%                            Passed to apply_mask. Default [] (no mask).
%     'parcellation'         atlas object. If supplied, returns an array of rsm
%                            objects, one per parcel. Default [] (no parcellation).
%     'smooth_mm'            Pre-smoothing FWHM in mm. Default [] (no smoothing).
%
%   Similarity-matrix options:
%     'treat_zero_as_data'   logical. Passed through to canlab_compute_similarity_matrix.
%                            Default false.
%
%   Misc:
%     'use_parallel'         logical. Default false.
%     'verbose'              logical. Default true.
%
% Output
% ------
%   R   rsm object, or [nParcels x 1] array of rsm objects when 'parcellation'
%       is set.

% =========================================================================
% Parse inputs
% =========================================================================
p = inputParser;
p.KeepUnmatched = false;
p.addParameter('group_by',            '',                @(x) ischar(x) || isstring(x) || iscellstr(x) || (iscell(x) && all(cellfun(@(c) ischar(c) || isstring(c), x))));
p.addParameter('condition_collapse',  'mean',            @(x) ischar(x) || isstring(x));
p.addParameter('level',               'subject',         @(x) ischar(x) || isstring(x));
p.addParameter('subject_var',         'subject_id',      @(x) ischar(x) || isstring(x));
p.addParameter('session_var',         'session_number',  @(x) ischar(x) || isstring(x));
p.addParameter('run_var',             'run_number',      @(x) ischar(x) || isstring(x));
p.addParameter('metric',              'correlation',     @(x) ischar(x) || isstring(x));
p.addParameter('fold_var',            'session_number',  @(x) ischar(x) || isstring(x));
p.addParameter('cv_scheme',           'allpairs',        @(x) (ischar(x)||isstring(x)) && ismember(lower(char(x)),{'allpairs','loo'}));
p.addParameter('whiten',              'none',            @(x) ischar(x) || isstring(x));
p.addParameter('whiten_method',       'covdiag',         @(x) ischar(x) || isstring(x));
p.addParameter('diagonal_correction', 'none',            @(x) ischar(x) || isstring(x));
p.addParameter('diagonal_group_by',   '',                @(x) ischar(x) || isstring(x));
p.addParameter('mask',                [],                @(x) true);
p.addParameter('parcellation',        [],                @(x) true);
p.addParameter('smooth_mm',           [],                @(x) isempty(x) || isnumeric(x));
p.addParameter('treat_zero_as_data',  false,             @(x) islogical(x) || isnumeric(x));
p.addParameter('nan_policy',          'propagate',       @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'propagate','skip_replicate'}));
p.addParameter('auto_groupings',      true,              @(x) islogical(x) || isnumeric(x));
p.addParameter('use_parallel',        false,             @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose',             true,              @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

% Normalize string args (leave group_by alone if cellstr — handled below)
fn = fieldnames(opt);
for i = 1:numel(fn)
    if strcmp(fn{i}, 'group_by'), continue; end
    if (ischar(opt.(fn{i})) || isstring(opt.(fn{i}))) && ~isempty(opt.(fn{i}))
        opt.(fn{i}) = char(opt.(fn{i}));
    end
end
% group_by: normalize to cellstr (length 1 = simple, length >1 = composite)
if isstring(opt.group_by), opt.group_by = cellstr(opt.group_by); end
if ischar(opt.group_by) && ~isempty(opt.group_by), opt.group_by = {opt.group_by}; end
if isempty(opt.group_by), opt.group_by = {}; end

% =========================================================================
% Validate combinations
% =========================================================================
% Fold-based (cross-validated) metrics require a fold column.
if ismember(lower(opt.metric), {'crossnobis','cvcorr','cvspearman'}) && isempty(opt.fold_var)
    error('compute_rsm:noFoldVar', ...
        'metric=''%s'' requires ''fold_var'' to identify the cross-validation fold column.', opt.metric);
end
if strcmpi(opt.metric, 'crossnobis')
    if strcmpi(opt.whiten, 'within_subject')
        warning('compute_rsm:doubleWhitening', ...
            ['Both crossnobis and within_subject whitening requested; crossnobis already ', ...
             'includes session-difference whitening. Skipping within_subject whitening.']);
        opt.whiten = 'session_difference';
    elseif strcmpi(opt.whiten, 'none')
        % Crossnobis without whitening is uncommon; default to session_difference
        opt.whiten = 'session_difference';
    end
end

% =========================================================================
% Pre-processing
% =========================================================================
if ~isempty(opt.smooth_mm)
    if opt.verbose, fprintf('compute_rsm: smoothing at FWHM = %g mm\n', opt.smooth_mm); end
    dat = preprocess(dat, 'smooth', opt.smooth_mm);
end

if ~isempty(opt.mask)
    if opt.verbose, fprintf('compute_rsm: applying mask\n'); end
    dat = apply_mask(dat, opt.mask);
end

% =========================================================================
% Parcellation branch — returns an array of rsm objects
% =========================================================================
if ~isempty(opt.parcellation)
    R = compute_rsm_parcelwise(dat, opt);
    return
end

% =========================================================================
% Single-mask branch
% =========================================================================
R = compute_rsm_one_mask(dat, opt, '');

end


% =========================================================================
function R_array = compute_rsm_parcelwise(dat, opt)

atlas_obj = opt.parcellation;
n_parcels = num_regions(atlas_obj);
labels    = atlas_obj.labels;
if numel(labels) < n_parcels
    labels = [labels, arrayfun(@(i) sprintf('parcel_%d', i), (numel(labels)+1):n_parcels, 'UniformOutput', false)];
end

R_array = rsm.empty;
if opt.verbose, fprintf('compute_rsm: running parcelwise (%d parcels)\n', n_parcels); end

% Resample the atlas into the DATA space ONCE (handles atlas/data space
% mismatch and avoids a per-parcel apply_mask, which is both slow and
% fragile). codes(v) gives the integer parcel index for data voxel v.
if compare_space(dat, atlas_obj) ~= 0
    if opt.verbose, fprintf('  resampling atlas to data space...\n'); end
    [~, atlas_rs] = evalc('resample_space(atlas_obj, dat)');
else
    atlas_rs = atlas_obj;
end
codes = round(double(atlas_rs.dat(:, 1)));

% Align parcel codes to dat.dat rows (drop removed voxels if present)
if ~isempty(dat.removed_voxels) && numel(dat.removed_voxels) == numel(codes)
    codes = codes(~dat.removed_voxels);
end
if numel(codes) ~= size(dat.dat, 1)
    % Last-resort fallback: per-parcel apply_mask (original slow path)
    if opt.verbose, warning('compute_rsm:parcelCodeMismatch', ...
        'Atlas codes (%d) do not align with data voxels (%d); using per-parcel apply_mask.', ...
        numel(codes), size(dat.dat,1)); end
    R_array = compute_rsm_parcelwise_applymask(dat, atlas_obj, labels, n_parcels, opt);
    return
end

dat_mat = double(dat.dat);   % [voxels x images]
for i = 1:n_parcels
    if opt.verbose && mod(i, 50) == 0, fprintf('  parcel %d/%d\n', i, n_parcels); end
    vox = (codes == i);
    if nnz(vox) < 2
        R_array(end+1) = rsm(); %#ok<AGROW>
        continue
    end
    dat_i = dat;
    dat_i.dat = dat_mat(vox, :);
    try
        R_array(end+1) = compute_rsm_one_mask(dat_i, opt, labels{i}); %#ok<AGROW>
    catch ME
        warning('compute_rsm:parcelFailed', 'Parcel %s failed: %s', labels{i}, ME.message);
        R_array(end+1) = rsm(); %#ok<AGROW>
    end
end

end


% =========================================================================
function R_array = compute_rsm_parcelwise_applymask(dat, atlas_obj, labels, n_parcels, opt)
% Fallback: per-parcel apply_mask (used only when integer-code alignment fails).
R_array = rsm.empty;
for i = 1:n_parcels
    try
        parcel_mask = select_atlas_subset(atlas_obj, labels(i));
        dat_i = apply_mask(dat, parcel_mask);
        R_array(end+1) = compute_rsm_one_mask(dat_i, opt, labels{i}); %#ok<AGROW>
    catch ME
        warning('compute_rsm:parcelFailed', 'Parcel %s failed: %s', labels{i}, ME.message);
        R_array(end+1) = rsm(); %#ok<AGROW>
    end
end
end


% =========================================================================
function R = compute_rsm_one_mask(dat, opt, parcel_name)
% Build a single rsm (possibly 3D across the replicate axis) for one mask.

mt = dat.metadata_table;
if isempty(mt) && ~strcmpi(opt.level, 'image')
    error('compute_rsm:noMetadata', ...
        'dat.metadata_table is empty but level=''%s'' was requested. Add metadata or use level=''image''.', opt.level);
end

% Determine replicate axis groupings
[replicate_groups, replicate_table] = determine_replicate_axis(mt, opt);
N = numel(replicate_groups);

% Determine condition (row/col) axis
[group_idx, group_labels, group_meta] = determine_group_axis(mt, opt);
k = numel(group_labels);

% Build per-replicate RSMs
dat_mat = double(dat.dat);   % [voxels x n_images]
out_3d  = nan(k, k, N);
whitened_info = struct('level', opt.whiten, 'method', opt.whiten_method, 'shrinkage', []);

% Quality-control accumulator: detect missing-condition cells per replicate.
% This is the upstream source of NaN in the output RSM stack.
qc = struct();
qc.nan_policy           = opt.nan_policy;
qc.n_replicates         = N;
qc.n_replicates_dropped = 0;
qc.replicate_qc         = repmat(struct( ...
    'replicate_index',         0, ...
    'replicate_label',         '', ...
    'n_missing_conditions',    0, ...
    'missing_condition_idx',   [], ...
    'missing_condition_labels', {{}}, ...
    'dropped',                 false), N, 1);

for n = 1:N
    rep_rows = replicate_groups{n};   % image indices belonging to this replicate
    if isempty(rep_rows), continue; end

    % Thread cross-validation fold_idx through opt (per-replicate).
    % Applies to all fold-based metrics: crossnobis + cvcorr/cvspearman.
    opt_n = opt;
    if ismember(lower(opt.metric), {'crossnobis','cvcorr','cvspearman'})
        if strcmpi(opt.fold_var, 'occurrence')
            % Auto-build folds by within-condition occurrence rank. Each
            % condition's 1st image -> fold 1, 2nd -> fold 2, etc. This aligns
            % every condition into the same fold set (the crossnobis
            % requirement) even when run/session numbering does not. Ordered
            % by run_number when present, else by row order.
            local_group = group_idx(rep_rows);
            if ismember(opt.run_var, mt.Properties.VariableNames)
                order_key = double(mt.(opt.run_var)(rep_rows));
            else
                order_key = (1:numel(rep_rows))';
            end
            fold_vals = zeros(numel(rep_rows), 1);
            ug = unique(local_group);
            for gg = 1:numel(ug)
                idx_g = find(local_group == ug(gg));
                [~, ord] = sort(order_key(idx_g));
                fold_vals(idx_g(ord)) = 1:numel(idx_g);
            end
            opt_n.fold_idx = fold_vals(:);
        else
            require_col(mt, opt.fold_var);
            fold_vals = mt.(opt.fold_var)(rep_rows);
            if iscategorical(fold_vals), fold_vals = double(fold_vals);
            elseif iscell(fold_vals) || isstring(fold_vals), [~, ~, fold_vals] = unique(fold_vals, 'stable');
            else, fold_vals = double(fold_vals);
            end
            opt_n.fold_idx = fold_vals(:);
        end
        if numel(unique(opt_n.fold_idx)) < 2
            warning('compute_rsm:tooFewFolds', ...
                'Replicate %d has only %d unique fold(s); %s returns NaN slice.', ...
                n, numel(unique(opt_n.fold_idx)), lower(opt.metric));
            out_3d(:, :, n) = NaN;
            continue
        end
    end

    [M, n_per_condition] = compute_one_rsm( ...
        dat_mat(:, rep_rows), ...
        group_idx(rep_rows), ...
        group_meta_subset(mt, rep_rows), ...
        k, opt_n);

    % Record per-replicate qc: which conditions had zero images
    qc.replicate_qc(n).replicate_index = n;
    if ~isempty(replicate_table)
        qc.replicate_qc(n).replicate_label = build_replicate_label(replicate_table, n);
    end
    if any(~isnan(n_per_condition))
        missing = find(n_per_condition == 0);
        qc.replicate_qc(n).n_missing_conditions   = numel(missing);
        qc.replicate_qc(n).missing_condition_idx  = missing;
        if ~isempty(group_labels)
            qc.replicate_qc(n).missing_condition_labels = group_labels(missing);
        end

        if ~isempty(missing) && strcmpi(opt.nan_policy, 'skip_replicate')
            qc.replicate_qc(n).dropped = true;
            qc.n_replicates_dropped = qc.n_replicates_dropped + 1;
            out_3d(:, :, n) = NaN;
            continue   % skip diagonal correction etc. for dropped slices
        end
    end

    % Diagonal correction
    if ~strcmpi(opt.diagonal_correction, 'none')
        diag_group_col = opt.diagonal_group_by;
        if isempty(diag_group_col), diag_group_col = opt.group_by; end
        if isempty(diag_group_col)
            warning('compute_rsm:noDiagGroup', 'diagonal_correction requested but no diagonal_group_by or group_by; skipping.');
        else
            % Build per-row diag-group at the k-condition level by mapping each
            % row to an integer ID based on its diag_group_col value. Resolve
            % the row's diag_group value from the first image of that condition,
            % then run a single findgroups() across ALL k rows so the resulting
            % integer IDs are consistent (rows with the same diag_group value
            % get the same integer).
            row_diag_vals = cell(k, 1);
            for g = 1:k
                first_idx = find(group_idx(rep_rows) == g, 1, 'first');
                if isempty(first_idx), row_diag_vals{g} = ''; continue; end
                abs_idx = rep_rows(first_idx);
                v = mt.(diag_group_col)(abs_idx);
                if iscell(v), v = v{1}; end
                if iscategorical(v) || isstring(v), v = char(v); end
                if isnumeric(v), v = num2str(v); end
                row_diag_vals{g} = char(v);
            end
            % Map distinct strings to integer group IDs
            [~, ~, row_diag_group] = unique(row_diag_vals, 'stable');
            if any(row_diag_group > 0)
                % The 'across_session_mean' / 'same_group_offdiag_mean' path
                % does not need per-row fold/session info -- pass NaN.
                M = apply_diagonal_correction(M, row_diag_group, [], opt.diagonal_correction);
            end
        end
    end

    out_3d(:, :, n) = M;
end

% qc summary across replicates
n_with_missing = sum(arrayfun(@(s) s.n_missing_conditions > 0, qc.replicate_qc));
qc.n_replicates_with_missing = n_with_missing;
if n_with_missing > 0 && opt.verbose
    if strcmpi(opt.nan_policy, 'skip_replicate')
        warning('compute_rsm:missingCells', ...
            ['%d of %d replicates had missing condition(s); %d slice(s) dropped (nan_policy=''skip_replicate''). ', ...
             'Inspect R.additional_info.qc.replicate_qc for per-replicate detail.'], ...
            n_with_missing, N, qc.n_replicates_dropped);
    else
        warning('compute_rsm:missingCells', ...
            ['%d of %d replicates had missing condition(s) -- those slices have NaN rows/cols ', ...
             '(nan_policy=''propagate''). Inspect R.additional_info.qc.replicate_qc, or re-run ', ...
             'with ''nan_policy'',''skip_replicate'' to drop bad slices.'], n_with_missing, N);
    end
end

% Within-subject and across-subject whitening operate on the stacked RSMs
if strcmpi(opt.whiten, 'within_subject')
    out_3d = whiten_rsm_stack(out_3d, opt.whiten_method);
    whitened_info.level = 'within_subject';
end
if strcmpi(opt.whiten, 'across_subject')
    out_3d = whiten_rsm_stack(out_3d, opt.whiten_method);
    whitened_info.level = 'across_subject';
end

% Determine is_dissimilarity from metric
is_dissim = ismember(lower(opt.metric), {'euclidean','seuclidean','mahalanobis','crossnobis'});

% Source string
if isempty(parcel_name)
    source = 'fmri_data';
else
    source = ['parcel:' parcel_name];
end

% Auto-build groupings from each group_by column (e.g., group_by={condition,
% bodysite} -> groupings.hot, groupings.warm, groupings.imagine,
% groupings.leftface, etc.) so downstream methods (cells, contrast,
% ttest_contrasts, reliability_by_grouping) work without manual setup.
auto_groupings = struct();
if logical(opt.auto_groupings) && ~isempty(opt.group_by) && ~isempty(group_meta) && k > 0
    auto_groupings = build_auto_groupings(group_meta, opt.group_by);
end

R = rsm(out_3d, ...
    'is_dissimilarity', is_dissim, ...
    'metric',           opt.metric, ...
    'labels',           group_labels, ...
    'metadata_table',   group_meta, ...
    'groupings',        auto_groupings, ...
    'level',            opt.level, ...
    'replicate_table',  replicate_table, ...
    'whitened',         whitened_info, ...
    'source',           source, ...
    'additional_info',  struct('qc', qc));

end


function groupings = build_auto_groupings(group_meta, group_by_cols)
% Build a bare-name struct of groupings from each group_by column. For
% group_meta with columns {condition, bodysite}, returns groupings like
% groupings.hot = [1:8], groupings.leftface = [1 9 17], etc.
%
% On name collisions across columns (rare), the later value wins; a single
% summary warning is emitted listing the conflicts.

groupings = struct();
collisions = {};
for c = 1:numel(group_by_cols)
    col = group_by_cols{c};
    if ~ismember(col, group_meta.Properties.VariableNames), continue; end
    v = group_meta.(col);
    % Get unique values preserving order
    if iscell(v) || iscategorical(v) || isstring(v)
        v_str = cellstr(string(v));
    else
        v_str = cellstr(string(v));
    end
    [unique_vals, ~, ~] = unique(v_str, 'stable');
    for j = 1:numel(unique_vals)
        name = sanitize_grouping_name(unique_vals{j});
        if isempty(name), continue; end
        idx = find(strcmp(v_str, unique_vals{j}));
        if isfield(groupings, name)
            collisions{end+1} = name; %#ok<AGROW>
        end
        groupings.(name) = idx(:)';
    end
end
if ~isempty(collisions)
    warning('compute_rsm:groupingCollision', ...
        ['Auto-groupings had name collisions on: %s. Later column values overwrote earlier. ', ...
         'Set ''auto_groupings'',false and attach R.groupings manually if you need both.'], ...
        strjoin(unique(collisions), ', '));
end
end


function s = sanitize_grouping_name(name)
% Strip leading/trailing whitespace, replace non-word characters with
% underscores, ensure leading letter so it's a valid struct field.
s = char(name);
s = strtrim(s);
s = regexprep(s, '[^A-Za-z0-9_]', '_');
if isempty(s), return; end
if ~isletter(s(1)), s = ['x' s]; end
end


function lbl = build_replicate_label(replicate_table, n)
% Build a short string describing replicate n from the table columns.
cols = replicate_table.Properties.VariableNames;
parts = cell(numel(cols), 1);
for i = 1:numel(cols)
    v = replicate_table.(cols{i})(n);
    parts{i} = sprintf('%s=%s', cols{i}, char(string(v)));
end
lbl = strjoin(parts, '/');
end


% =========================================================================
function [groups, replicate_table] = determine_replicate_axis(mt, opt)

switch lower(opt.level)
    case 'collapsed'
        groups = {(1:height(mt))'};
        replicate_table = table({'all'}, 'VariableNames', {'aggregation'});
    case 'subject'
        require_col(mt, opt.subject_var);
        [G, ID] = findgroups(mt.(opt.subject_var));
        groups = arrayfun(@(g) find(G == g), 1:max(G), 'UniformOutput', false)';
        replicate_table = table(ID, 'VariableNames', {opt.subject_var});
    case 'session'
        require_col(mt, opt.subject_var);
        require_col(mt, opt.session_var);
        [G, S, Sess] = findgroups(mt.(opt.subject_var), mt.(opt.session_var));
        groups = arrayfun(@(g) find(G == g), 1:max(G), 'UniformOutput', false)';
        replicate_table = table(S, Sess, 'VariableNames', {opt.subject_var, opt.session_var});
    case 'run'
        require_col(mt, opt.subject_var);
        require_col(mt, opt.session_var);
        require_col(mt, opt.run_var);
        [G, S, Sess, Run] = findgroups(mt.(opt.subject_var), mt.(opt.session_var), mt.(opt.run_var));
        groups = arrayfun(@(g) find(G == g), 1:max(G), 'UniformOutput', false)';
        replicate_table = table(S, Sess, Run, 'VariableNames', {opt.subject_var, opt.session_var, opt.run_var});
    case 'image'
        % One slice per subject, but no condition collapsing — image-level RSM per subject
        require_col(mt, opt.subject_var);
        [G, ID] = findgroups(mt.(opt.subject_var));
        groups = arrayfun(@(g) find(G == g), 1:max(G), 'UniformOutput', false)';
        replicate_table = table(ID, 'VariableNames', {opt.subject_var});
    otherwise
        error('compute_rsm:badLevel', 'level must be one of {subject, session, run, collapsed, image}; got %s.', opt.level);
end

end


% =========================================================================
function [group_idx, group_labels, group_meta] = determine_group_axis(mt, opt)
% Returns:
%   group_idx     [n_images x 1] integer label (1..k) per image
%   group_labels  {k x 1} cellstr names for each row/col
%   group_meta    [k x 1] table of metadata for each row/col
%
% opt.group_by is a cellstr (length 1 = simple grouping, length >1 = composite).

if isempty(opt.group_by) || strcmpi(opt.level, 'image')
    % Image-level: each row/col is its own image
    n = height(mt);
    group_idx = (1:n)';
    if n > 0
        group_labels = cellfun(@(i) sprintf('img_%d', i), num2cell(1:n)', 'UniformOutput', false);
        group_meta = mt;
    else
        group_labels = {};
        group_meta = mt;
    end
    return
end

cols = opt.group_by;
for i = 1:numel(cols), require_col(mt, cols{i}); end

if numel(cols) == 1
    v = mt.(cols{1});
    [group_idx, G_vals] = findgroups(v);
    group_labels = cellstr(string(G_vals));
else
    % Composite grouping: stringify each column and concatenate with '_'
    parts = cell(height(mt), numel(cols));
    for i = 1:numel(cols)
        v = mt.(cols{i});
        if iscell(v)
            parts(:, i) = cellfun(@(x) char(string(x)), v, 'UniformOutput', false);
        else
            parts(:, i) = cellstr(string(v));
        end
    end
    composite_key = strings(height(mt), 1);
    for r = 1:height(mt)
        composite_key(r) = strjoin(parts(r, :), '_');
    end
    [group_idx, G_vals] = findgroups(composite_key);
    group_labels = cellstr(string(G_vals));
end

% Per-condition metadata: take first row per group
group_meta = table.empty;
if ~isempty(mt)
    k = numel(group_labels);
    rows = cell(k, 1);
    for g = 1:k
        idxs = find(group_idx == g);
        rows{g} = mt(idxs(1), :);
    end
    group_meta = vertcat(rows{:});
end

end


% =========================================================================
function [M, n_per_condition] = compute_one_rsm(X, group_idx, ~, k, opt)
% X         [voxels x n_images] data
% group_idx [n_images x 1] integer labels (1..k)
% k         total number of groups (rows/cols in output)
% opt       option struct
%
% Returns:
%   M               [k x k] similarity / dissimilarity
%   n_per_condition [k x 1] image counts per condition. NaN-filled for the
%                   crossnobis path (which uses fold-level aggregation).

n_per_condition = nan(k, 1);

% --- Crossnobis path ---
if strcmpi(opt.metric, 'crossnobis')
    M = compute_crossnobis(X, group_idx, opt);
    return
end

% --- Cross-validated (cross-fold) correlation path ---
if ismember(lower(opt.metric), {'cvcorr','cvspearman'})
    M = compute_cvcorr(X, group_idx, opt);
    return
end

% --- All other metrics path ---
% Collapse to one pattern per condition
[P, n_per_condition] = aggregate_patterns(X, group_idx, k, opt.condition_collapse);   % [k x voxels]

switch lower(opt.metric)
    case 'correlation'
        % Use pairwise-missing-aware similarity primitive
        try
            M = canlab_compute_similarity_matrix(P', ...
                'similarity_metric', 'correlation', ...
                'treat_zero_as_data', logical(opt.treat_zero_as_data), ...
                'verbose', false);
        catch
            M = corr(P', 'Type', 'Pearson', 'rows', 'pairwise');
        end
    case 'spearman'
        % Routed through the canlab primitive (Spearman support added 2026-05).
        try
            M = canlab_compute_similarity_matrix(P', ...
                'similarity_metric', 'spearman', ...
                'treat_zero_as_data', logical(opt.treat_zero_as_data), ...
                'verbose', false);
        catch
            M = corr(P', 'Type', 'Spearman', 'rows', 'pairwise');
        end
    case 'cosine'
        try
            M = canlab_compute_similarity_matrix(P', ...
                'similarity_metric', 'cosine_similarity', ...
                'treat_zero_as_data', logical(opt.treat_zero_as_data), ...
                'verbose', false);
        catch
            d = pdist(P, 'cosine');
            M = 1 - squareform(d);
            M(1:k+1:end) = 1;
        end
    case 'euclidean'
        d = pdist(P, 'euclidean');
        M = squareform(d);
    case 'seuclidean'
        d = pdist(P, 'seuclidean');
        M = squareform(d);
    case 'mahalanobis'
        % Regularized covariance from the rows of P
        try
            sigma = ledoit_wolf_shrinkage(P);
            d = pdist(P, 'mahalanobis', sigma);
        catch
            d = pdist(P, 'mahalanobis');
        end
        M = squareform(d);
    otherwise
        error('compute_rsm:badMetric', 'Unknown metric: %s', opt.metric);
end

end


% =========================================================================
function M = compute_crossnobis(X, group_idx, opt)
% X         [voxels x n_images]
% group_idx [n_images x 1] integer condition labels (1..k)

% Need fold labels for this slice. Look up from the parent function's
% metadata via a side channel: we re-derive from group_idx + opt.fold_var
% by relying on the fact that compute_rsm_one_mask passes us the per-image
% group_idx but NOT the fold_idx. We need fold_idx here.
%
% Solution: stash the fold_idx in opt at the caller before invoking.
% But for clarity, we instead require the caller to put fold_idx in opt.
% Since opt comes from compute_rsm's parse — wire fold_idx through there.

if ~isfield(opt, 'fold_idx') || isempty(opt.fold_idx)
    error('compute_rsm:crossnobisNoFold', ...
        'Internal error: crossnobis path requires opt.fold_idx to be set by compute_rsm_one_mask.');
end

fold_idx = opt.fold_idx;   % [n_images x 1]

% Build per-fold pattern matrices [k x voxels]
Xfolds = build_fold_pattern_matrices(X, group_idx, fold_idx);

% Whiten via session-difference (the standard crossnobis pre-whitening)
if any(strcmpi(opt.whiten, {'session_difference','within_subject','across_subject'}))
    Xfolds = whiten_session_difference(Xfolds);
end

% Mean-center across conditions within each fold
for f = 1:numel(Xfolds)
    Xfolds{f} = Xfolds{f} - mean(Xfolds{f}, 1, 'omitnan');
end

% Cross-validated dissimilarity: mean over all distinct fold pairs of
% (X_i[a] - X_i[b]) .* (X_j[a] - X_j[b])
[k, ~] = size(Xfolds{1});
M = zeros(k, k);
n_folds = numel(Xfolds);
fold_pairs = nchoosek(1:n_folds, 2);
n_pairs = size(fold_pairs, 1);

for a = 1:k
    for b = (a+1):k
        d_accum = 0;
        n_used = 0;
        for p = 1:n_pairs
            i = fold_pairs(p, 1); j = fold_pairs(p, 2);
            di = Xfolds{i}(a, :) - Xfolds{i}(b, :);
            dj = Xfolds{j}(a, :) - Xfolds{j}(b, :);
            mask = isfinite(di) & isfinite(dj);
            if any(mask)
                d_accum = d_accum + mean(di(mask) .* dj(mask));
                n_used = n_used + 1;
            end
        end
        if n_used > 0
            M(a, b) = d_accum / n_used;
            M(b, a) = M(a, b);
        end
    end
end

end


% =========================================================================
function M = compute_cvcorr(X, group_idx, opt)
% Cross-validated (cross-fold) correlation SIMILARITY matrix.
%
% The two condition patterns being correlated ALWAYS come from different
% folds, so any within-fold shared structure (e.g. the shared-run noise +
% run/bodysite mean when Hot/Warm/Imagine co-occur in one run) does NOT
% contribute to the correlation. This reproduces the Sun et al. (2026)
% BodyMap RSM construction (per-run RSMs, exclude within-session/run
% correlations, average the rest) as a single metric. Unlike crossnobis it
% stays in correlation (similarity) space; the diagonal M(i,i) is the
% cross-fold reliability of condition i (NOT 1).
%
% cv_scheme:
%   'allpairs' (default) -- M(i,j) = mean over ordered fold pairs (a,b), a~=b,
%       of corr(Xfolds{a}(i,:), Xfolds{b}(j,:)). Most conservative: every
%       correlation is between two SINGLE-fold patterns (noisiest, lowest
%       power, fully symmetric use of the data).
%   'loo' -- leave-one-fold-out. For each held-out fold f,
%       M(i,j) += corr(Xfolds{f}(i,:), meanOfRest_f(j,:)), where meanOfRest_f
%       is the mean over the OTHER folds. The training side is averaged over
%       n_folds-1 folds => higher SNR, less attenuation, MORE POWER, while
%       still never sharing a fold. Closer to the paper's whole-brain
%       leave-session-out masking. Result symmetrized over held-out/train roles.
%
% metric='cvcorr' -> Pearson; metric='cvspearman' -> Spearman.

if ~isfield(opt, 'fold_idx') || isempty(opt.fold_idx)
    error('compute_rsm:cvcorrNoFold', ...
        'Internal error: cvcorr path requires opt.fold_idx to be set by compute_rsm.');
end

Xfolds  = build_fold_pattern_matrices(X, group_idx, opt.fold_idx);
n_folds = numel(Xfolds);
k       = size(Xfolds{1}, 1);

if strcmpi(opt.metric, 'cvspearman'), ctype = 'Spearman'; else, ctype = 'Pearson'; end
scheme = 'allpairs';
if isfield(opt, 'cv_scheme') && ~isempty(opt.cv_scheme), scheme = lower(char(opt.cv_scheme)); end

M_accum = zeros(k, k);
M_count = zeros(k, k);

switch scheme
    case 'allpairs'
        for a = 1:n_folds
            for b = 1:n_folds
                if a == b, continue; end
                % Rab(i,j) = corr(condition i in fold a, condition j in fold b)
                Rab = corr(Xfolds{a}', Xfolds{b}', 'type', ctype, 'rows', 'pairwise');
                good = isfinite(Rab);
                M_accum(good) = M_accum(good) + Rab(good);
                M_count(good) = M_count(good) + 1;
            end
        end

    case 'loo'
        % Per-condition running sum + present-count across folds so the
        % mean-of-rest can be formed by subtracting the held-out fold. Handles
        % conditions missing from some folds (partial sessions).
        n_vox = size(Xfolds{1}, 2);
        Xsum  = zeros(k, n_vox);
        pres  = zeros(k, 1);
        for f = 1:n_folds
            Xf = Xfolds{f};
            pf = all(isfinite(Xf), 2);   % [k x 1] conditions present in fold f
            Xf0 = Xf; Xf0(~pf, :) = 0;
            Xsum = Xsum + Xf0;
            pres = pres + pf;
        end
        for f = 1:n_folds
            Xf = Xfolds{f};
            pf = all(isfinite(Xf), 2);
            Xf0 = Xf; Xf0(~pf, :) = 0;
            rest_cnt = pres - double(pf);          % folds contributing to rest, per condition
            Rest = nan(k, n_vox);
            v = rest_cnt > 0;
            Rest(v, :) = (Xsum(v, :) - Xf0(v, :)) ./ rest_cnt(v);
            % Rf(i,j) = corr(held-out fold cond i, mean-of-rest cond j)
            Rf = corr(Xf', Rest', 'type', ctype, 'rows', 'pairwise');
            good = isfinite(Rf);
            M_accum(good) = M_accum(good) + Rf(good);
            M_count(good) = M_count(good) + 1;
        end

    otherwise
        error('compute_rsm:badCvScheme', ...
            'cv_scheme must be ''allpairs'' or ''loo''; got %s.', scheme);
end

M = M_accum ./ max(M_count, 1);
M(M_count == 0) = NaN;
M = (M + M') / 2;         % symmetric in expectation / over held-out roles

end


% =========================================================================
function [P, n_per_condition] = aggregate_patterns(X, group_idx, k, mode)
% Average (or otherwise aggregate) X within each group, returning [k x voxels].
%
% Second output n_per_condition is a [k x 1] integer vector of how many
% images contributed to each group. Used by compute_rsm_one_mask to detect
% missing-condition replicates for the nan_policy logic.

n_vox = size(X, 1);
P = nan(k, n_vox);
n_per_condition = zeros(k, 1);

switch lower(mode)
    case 'mean'
        for g = 1:k
            rows = group_idx == g;
            n_per_condition(g) = nnz(rows);
            if any(rows)
                P(g, :) = mean(X(:, rows), 2, 'omitnan')';
            end
        end
    case 'concat'
        error('compute_rsm:concatNotImplemented', ...
            '''concat'' aggregation is reserved for image-level RSMs (level=''image''); not implemented here.');
    case 'none'
        if any(accumarray(group_idx(:), 1) > 1)
            error('compute_rsm:noneRequiresSingleton', ...
                'condition_collapse=''none'' requires exactly one image per condition.');
        end
        for g = 1:k
            row = find(group_idx == g, 1);
            P(g, :) = X(:, row)';
        end
    otherwise
        error('compute_rsm:badCollapse', 'Unknown condition_collapse: %s', mode);
end

end


% =========================================================================
function out_3d = whiten_rsm_stack(out_3d, method)
% Vectorize upper triangle of each k x k slice -> [N x p] matrix,
% whiten across rows, re-square.

[k, ~, N] = size(out_3d);
mask = triu(true(k), 1);
p = nnz(mask);
X = zeros(N, p);
for n = 1:N
    slice = out_3d(:, :, n);
    X(n, :) = slice(mask)';
end

Xw = whiten_within_subject(X, method);

for n = 1:N
    slice = zeros(k);
    slice(mask) = Xw(n, :);
    slice = slice + slice';
    slice(1:k+1:end) = 1;   % diagonal heuristic: 1 for RSM, ignored for RDM display
    out_3d(:, :, n) = slice;
end

end


% =========================================================================
function require_col(mt, col)
if ~ismember(col, mt.Properties.VariableNames)
    error('compute_rsm:missingColumn', ...
        'dat.metadata_table is missing required column ''%s''.', col);
end
end


function v = coerce_to_int(x)
if iscategorical(x), v = double(x);
elseif iscell(x), [~, ~, v] = unique(x, 'stable');
elseif isstring(x), [~, ~, v] = unique(x, 'stable');
else, v = double(x);
end
v = v(1);
end


function meta_sub = group_meta_subset(mt, rep_rows) %#ok<INUSD>
% Currently unused; reserved for future condition_collapse='none' / merging paths
meta_sub = [];
end
