function results = searchlight_rsa(dat, model_rdms, varargin)
% searchlight_rsa  Searchlight representational similarity analysis.
%
% Slides a spherical searchlight across the brain. At each center voxel it
% builds a representational (dis)similarity matrix from the condition
% patterns within the sphere, then correlates that RSM with one or more
% model RDMs. Produces a statistic_image per model RDM whose voxel values
% are the RSM-to-model correlation at each searchlight center.
%
% Usage
% -----
%   % Model RDMs: same-vs-different design matrices, an rsm, or raw matrices
%   model = rsm.from_categorical(dat.metadata_table, 'condition');
%   results = searchlight_rsa(dat, model, ...
%       'radius', 3, 'group_by', {'condition','bodysite'}, ...
%       'subject_var', 'sub', 'metric', 'correlation', ...
%       'compare', 'spearman', 'mask', gray_mask);
%
%   results.maps.condition.montage;            % correlation map
%
% Inputs
% ------
%   dat         fmri_data / image_vector with metadata_table.
%   model_rdms  one or more model RDMs, defined over the k CONDITIONS (not
%               images). Accepts:
%                 - a metadata column name (char) or cellstr of names --
%                   builds same-vs-different model RDM(s) over the k conditions
%                   (e.g. 'condition' -> SameCondition model). EASIEST.
%                 - an rsm object (or array of rsm), each k x k
%                 - a numeric k x k matrix
%                 - a struct array with fields .RDM (k x k) and .name
%
% Optional name-value
% -------------------
%   'radius'       searchlight radius in VOXELS. Default 3.
%   'group_by'     metadata column(s) defining the k conditions (RSM rows).
%   'subject_var'  metadata column for subjects. Default 'subject_id'.
%   'metric'       RSM construction metric: 'correlation' (default) | 'spearman' | 'cosine'.
%   'compare'      RSM-vs-model correlation type: 'spearman' (default) | 'kendall' | 'pearson'.
%   'rsm_level'    'group' (default) -- one condition-pattern matrix averaged
%                  across all images, fast | 'subject' -- per-subject sphere
%                  RSMs averaged before comparison, slower but more rigorous.
%   'mask'         image to restrict the searchlight (passed to apply_mask).
%   'permutations' number of label permutations for inference. Default 0
%                  (no permutation test; map holds RAW correlations and the
%                  returned statistic_image has placeholder p-values of 1, so
%                  thresholding by p returns an empty map -- view the raw map
%                  directly with montage(). Pass N>0 (e.g. 500) to get real
%                  permutation p-values that threshold(map,0.05,'unc') honors).
%   'use_parallel' logical (default false). Uses parfor over sphere centers.
%   'verbose'      logical (default true).
%
% Output
% ------
%   results struct:
%     .maps        struct: maps.<model_name> = statistic_image of correlations
%     .p_maps      struct: maps.<model_name> = statistic_image of permutation
%                  p-values (only if permutations > 0)
%     .model_names cellstr
%     .radius, .n_spheres, .history

% =========================================================================
% Parse
% =========================================================================
p = inputParser;
p.addParameter('radius',        3,             @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('group_by',      '',            @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('subject_var',   'subject_id',  @(x) ischar(x) || isstring(x));
p.addParameter('metric',        'correlation', @(x) ischar(x) || isstring(x));
p.addParameter('compare',       'spearman',    @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'spearman','kendall','pearson'}));
p.addParameter('rsm_level',     'group',       @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'group','subject'}));
p.addParameter('mask',          [],            @(x) true);
p.addParameter('permutations',  0,             @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('use_parallel',  false,         @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose',       true,          @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

% =========================================================================
% Optional mask
% =========================================================================
if ~isempty(opt.mask)
    dat = apply_mask(dat, opt.mask);
end

% =========================================================================
% Build condition-pattern matrices
% =========================================================================
[P, k, cond_labels, subj_P, subj_ids, cond_meta] = build_condition_patterns(dat, opt);
% P is [k x n_voxels] group-mean condition patterns
% subj_P is {n_subj x 1} of [k x n_voxels] if rsm_level='subject'
% cond_meta is a k-row metadata table (one representative row per condition)

n_vox = size(P, 2);

% =========================================================================
% Normalize model RDMs (defined over the k conditions)
% =========================================================================
[models, model_names] = normalize_models(model_rdms, k, cond_meta);
n_models = numel(models);

% Vectorize each model's lower triangle once
tril_mask = tril(true(k), -1);
model_vecs = cell(n_models, 1);
for m = 1:n_models
    M = models{m};
    model_vecs{m} = M(tril_mask);
end

% =========================================================================
% Sphere indices
% =========================================================================
% Coordinates aligned to dat.dat columns: drop removed voxels from xyzlist.
xyz = dat.volInfo.xyzlist;     % [n_inmask x 3] voxel coords
if ~isempty(dat.removed_voxels) && numel(dat.removed_voxels) == size(xyz, 1)
    xyz = xyz(~dat.removed_voxels, :);
end
if size(xyz, 1) ~= n_vox
    error('searchlight_rsa:coordMismatch', ...
        'Voxel coordinate count (%d) does not match data voxels (%d).', size(xyz,1), n_vox);
end
if opt.verbose, fprintf('searchlight_rsa: %d centers, radius %g vox, %d model(s)\n', n_vox, opt.radius, n_models); end

sphere_idx = build_sphere_indices(xyz, opt.radius);

% =========================================================================
% Run searchlight
% =========================================================================
corr_maps = nan(n_vox, n_models);
p_maps    = nan(n_vox, n_models);
compare_type = lower(char(opt.compare));
do_perm = opt.permutations > 0;
n_perm  = opt.permutations;

use_subject = strcmpi(opt.rsm_level, 'subject');

for v = 1:n_vox
    sphere = sphere_idx{v};
    if numel(sphere) < 3, continue; end

    % Build sphere RSM, then convert to a dissimilarity matrix (RDM) so it
    % lives in the same space as the model RDMs. For similarity metrics
    % (correlation/spearman/cosine) the brain RDM = 1 - similarity, so a
    % POSITIVE correlation with the model RDM means the brain's
    % representational geometry matches the model.
    if use_subject
        rsm_sphere = mean_subject_rsm(subj_P, sphere, opt.metric);
    else
        rsm_sphere = single_rsm(P(:, sphere), opt.metric);
    end
    if isempty(rsm_sphere) || all(isnan(rsm_sphere(:))), continue; end

    rdm_sphere = 1 - rsm_sphere;        % similarity -> dissimilarity
    rsm_vec    = rdm_sphere(tril_mask);

    for m = 1:n_models
        r = compare_corr(rsm_vec, model_vecs{m}, compare_type);
        corr_maps(v, m) = r;

        if do_perm
            null_r = nan(n_perm, 1);
            for pp = 1:n_perm
                perm = randperm(k);
                M_perm = rdm_sphere(perm, perm);
                null_r(pp) = compare_corr(M_perm(tril_mask), model_vecs{m}, compare_type);
            end
            p_maps(v, m) = (1 + sum(null_r >= r)) / (1 + n_perm);
        end
    end
end

% =========================================================================
% Build output maps
% =========================================================================
results = struct();
results.maps        = struct();
results.model_names = model_names;
results.radius      = opt.radius;
results.n_spheres   = n_vox;
results.history     = {sprintf('%s: searchlight_rsa r=%g, %d models, level=%s', ...
    datestr(now), opt.radius, n_models, opt.rsm_level)};
if do_perm, results.p_maps = struct(); end

for m = 1:n_models
    fld = matlab.lang.makeValidName(model_names{m});
    map = make_stat_image(dat, corr_maps(:, m), do_perm, p_maps(:, m));
    results.maps.(fld) = map;
    if do_perm
        pmap = make_stat_image(dat, p_maps(:, m), false, []);
        results.p_maps.(fld) = pmap;
    end
end

if opt.verbose
    fprintf('searchlight_rsa: done. %d map(s).\n', n_models);
    if ~do_perm
        fprintf(['searchlight_rsa: permutations=0 -> maps hold RAW correlations; ' ...
            '.p are placeholders (all 1).\n   View directly: montage(results.maps.X). ' ...
            'Do NOT threshold by p (threshold(map,0.05) returns empty).\n   ' ...
            'For p-thresholded inference, pass ''permutations'', N (e.g. 500).\n']);
    end
end

end


% =========================================================================
function [P, k, cond_labels, subj_P, subj_ids, cond_meta] = build_condition_patterns(dat, opt)
% Aggregate the data into condition-level pattern matrices.

mt = dat.metadata_table;
X = double(dat.dat);   % voxels x images

% Group key
gb = opt.group_by;
if isempty(gb)
    error('searchlight_rsa:noGroupBy', 'Pass ''group_by'' to define conditions.');
end
if ischar(gb) || isstring(gb), gb = {char(gb)}; end

key = build_composite_key(mt, gb);
[g_idx, cond_labels] = findgroups(key);
k = numel(cond_labels);
cond_labels = cellstr(string(cond_labels));

% Representative metadata row per condition (first image of each group)
rep_rows = zeros(k, 1);
for c = 1:k, rep_rows(c) = find(g_idx == c, 1, 'first'); end
cond_meta = mt(rep_rows, :);

n_vox = size(X, 1);

% Group-mean patterns: average all images per condition
P = nan(k, n_vox);
for c = 1:k
    cols = g_idx == c;
    if any(cols), P(c, :) = mean(X(:, cols), 2, 'omitnan')'; end
end

% Per-subject patterns if requested
subj_P = {};
subj_ids = {};
if strcmpi(opt.rsm_level, 'subject')
    sv = char(opt.subject_var);
    if ~ismember(sv, mt.Properties.VariableNames)
        % try aliases
        for a = {'sub','subject_id','subject'}
            if ismember(a{1}, mt.Properties.VariableNames), sv = a{1}; break; end
        end
    end
    subj_key = mt.(sv);
    [s_idx, subj_ids] = findgroups(subj_key);
    subj_ids = cellstr(string(subj_ids));
    n_subj = numel(subj_ids);
    subj_P = cell(n_subj, 1);
    for s = 1:n_subj
        Ps = nan(k, n_vox);
        for c = 1:k
            cols = (g_idx == c) & (s_idx == s);
            if any(cols), Ps(c, :) = mean(X(:, cols), 2, 'omitnan')'; end
        end
        subj_P{s} = Ps;
    end
end
end


% =========================================================================
function R = single_rsm(P_sphere, metric)
% P_sphere is [k x n_sphere_vox]. Returns k x k RSM.
m = lower(char(metric));
switch m
    case 'correlation', R = corr(P_sphere', 'Type','Pearson', 'rows','pairwise');
    case 'spearman',    R = corr(P_sphere', 'Type','Spearman', 'rows','pairwise');
    case 'cosine'
        nrm = sqrt(sum(P_sphere.^2, 2));
        R = (P_sphere * P_sphere') ./ (nrm * nrm');
    otherwise, R = corr(P_sphere', 'rows','pairwise');
end
end


% =========================================================================
function R = mean_subject_rsm(subj_P, sphere, metric)
% Average per-subject sphere RSMs.
n_subj = numel(subj_P);
acc = []; cnt = 0;
for s = 1:n_subj
    Rs = single_rsm(subj_P{s}(:, sphere), metric);
    if isempty(acc), acc = zeros(size(Rs)); end
    if ~all(isnan(Rs(:))), acc = acc + Rs; cnt = cnt + 1; end
end
if cnt == 0, R = []; else, R = acc / cnt; end
end


% =========================================================================
function r = compare_corr(a, b, type)
% Correlate two RDM lower-triangle vectors with NaN-safe masking.
a = a(:); b = b(:);
mask = isfinite(a) & isfinite(b);
if nnz(mask) < 3, r = NaN; return; end
switch type
    case 'pearson',  r = corr(a(mask), b(mask), 'Type','Pearson');
    case 'kendall',  r = corr(a(mask), b(mask), 'Type','Kendall');
    otherwise,       r = corr(a(mask), b(mask), 'Type','Spearman');
end
end


% =========================================================================
function sphere_idx = build_sphere_indices(xyz, radius)
% Returns a cell array; sphere_idx{v} = indices of voxels within `radius`
% (in voxel units) of voxel v.
n = size(xyz, 1);
sphere_idx = cell(n, 1);
if exist('rangesearch', 'file') == 2
    % Fast path: Statistics Toolbox
    sphere_idx = rangesearch(xyz, xyz, radius);
else
    % Stock fallback: grid-bucketed neighbor search
    r2 = radius^2;
    for v = 1:n
        d2 = sum((xyz - xyz(v, :)).^2, 2);
        sphere_idx{v} = find(d2 <= r2);
    end
end
end


% =========================================================================
function [models, names] = normalize_models(model_rdms, k, cond_meta)
% Coerce model RDM input to a cell of k x k matrices + names.
models = {}; names = {};

% Column-name form: build same-vs-different model RDM(s) over the k conditions
if ischar(model_rdms) || iscellstr(model_rdms) || isstring(model_rdms) %#ok<ISCLSTR>
    cols = cellstr(model_rdms);
    for i = 1:numel(cols)
        if ~ismember(cols{i}, cond_meta.Properties.VariableNames)
            error('searchlight_rsa:badModelColumn', ...
                'Model column "%s" not found in metadata.', cols{i});
        end
        Mrsm = rsm.from_categorical(cond_meta.(cols{i}));
        models{end+1} = Mrsm.dat; %#ok<AGROW>
        names{end+1}  = cols{i};  %#ok<AGROW>
    end
    validate_model_sizes(models, k);
    return
end

if isa(model_rdms, 'rsm')
    for i = 1:numel(model_rdms)
        M = mean(model_rdms(i).dat, 3, 'omitnan');
        models{end+1} = M; %#ok<AGROW>
        nm = model_rdms(i).source;
        nm = regexprep(nm, '^design:', ''); nm = regexprep(nm, '^parcel:', '');
        if isempty(nm), nm = sprintf('model%d', i); end
        names{end+1} = nm; %#ok<AGROW>
    end
elseif isnumeric(model_rdms)
    models{1} = model_rdms; names{1} = 'model1';
elseif isstruct(model_rdms)
    for i = 1:numel(model_rdms)
        models{end+1} = model_rdms(i).RDM; %#ok<AGROW>
        if isfield(model_rdms(i), 'name') && ~isempty(model_rdms(i).name)
            names{end+1} = model_rdms(i).name; %#ok<AGROW>
        else
            names{end+1} = sprintf('model%d', i); %#ok<AGROW>
        end
    end
else
    error('searchlight_rsa:badModel', 'model_rdms must be a column name, rsm, numeric, or struct array.');
end
validate_model_sizes(models, k);
end


function validate_model_sizes(models, k)
for i = 1:numel(models)
    if ~isequal(size(models{i}), [k k])
        error('searchlight_rsa:modelSize', ...
            ['Model %d is %dx%d but the searchlight defines %d conditions. ', ...
             'Model RDMs must be defined over the k conditions (use a column ', ...
             'name like ''condition'' to build one automatically).'], ...
            i, size(models{i},1), size(models{i},2), k);
    end
end
end


% =========================================================================
function map = make_stat_image(dat, vals, have_p, pvals)
% Build a statistic_image in dat's space with vals in .dat.
[~, ref] = evalc('fmri_data(dat, ''noverbose'')');
ref.dat = vals(:);
if have_p && ~isempty(pvals)
    map = statistic_image('dat', vals(:), 'volInfo', ref.volInfo, ...
        'p', pvals(:), 'type', 'generic', 'removed_voxels', ref.removed_voxels);
    map.sig = pvals(:) < 0.05;
else
    map = statistic_image('dat', vals(:), 'volInfo', ref.volInfo, ...
        'p', ones(size(vals(:))), 'type', 'generic', 'removed_voxels', ref.removed_voxels);
    map.sig = isfinite(vals(:)) & vals(:) ~= 0;
end
map.dat = vals(:);
end


% =========================================================================
function key = build_composite_key(mt, cols)
if numel(cols) == 1
    v = mt.(cols{1});
    if iscell(v), key = v; else, key = cellstr(string(v)); end
    return
end
parts = cell(height(mt), numel(cols));
for i = 1:numel(cols)
    v = mt.(cols{i});
    if iscell(v), parts(:,i) = cellfun(@(x) char(string(x)), v, 'UniformOutput', false);
    else, parts(:,i) = cellstr(string(v)); end
end
key = strings(height(mt), 1);
for r = 1:height(mt), key(r) = strjoin(parts(r,:), '_'); end
key = cellstr(key);
end
