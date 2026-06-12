function results = rsa_parcelwise(dat, varargin)
% rsa_parcelwise  End-to-end parcelwise RSA emitting statistic_image brain maps.
%
% Builds a per-parcel RSM for every region in an atlas, runs RSA inference
% (declarative contrasts and/or a per-parcel LME), FDR-corrects across
% parcels, and projects the results back onto the brain as statistic_image
% objects so the standard .threshold / .montage / .table chain just works.
%
% Reproduces the Sun et al. workflow's analyze_roi_contrasts + region2fmri_data
% pipeline (Figs 7G, S8, S9) in a single call.
%
% Usage
% -----
%   atlas = load_atlas('canlab2024');
%
%   % (A) Contrast path: within/between condition contrasts per parcel
%   spec = {
%       'hot_vs_warm',   'hot',   'warm';
%       'HW_vs_HI',      {'hot','warm'}, {'hot','imagine'};
%   };
%   results = rsa_parcelwise(dat, 'atlas', atlas, ...
%       'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
%       'metric', 'spearman', 'contrasts', spec, ...
%       'correction', 'fdr', 'tail', 'right');
%
%   results.maps.hot_vs_warm.threshold(0.05, 'unc').montage;
%
%   % (B) LME path: per-parcel mixed-effects model
%   results = rsa_parcelwise(dat, 'atlas', atlas, ...
%       'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
%       'lme', 'Y ~ SameCondition + SameBodysite + (1|Sub)', ...
%       'predictors', {'condition','bodysite'});
%
%   results.maps.SameCondition.montage;
%
% Inputs
% ------
%   dat   fmri_data object with metadata_table.
%
% Required name-value
% -------------------
%   'atlas'   atlas object defining the parcels.
%
% Aggregation / metric (forwarded to compute_rsm)
%   'group_by', 'subject_var', 'session_var', 'metric', 'level',
%   'nan_policy', 'fold_var', 'whiten' -- see compute_rsm.
%
% Analysis (at least one of contrasts / lme required)
%   'contrasts'  ttest_contrasts spec (Nx3 cell). Runs per parcel.
%   'lme'        Wilkinson formula string for a per-parcel rsa_lme. Requires
%                'predictors' to be passed too.
%   'predictors' cellstr of metadata columns for the LME path.
%   'interactions' / 'three_way'  forwarded to rsa_lme.
%
% Inference options
%   'tail'        'both' (default) | 'right' | 'left'   (contrast path)
%   'correction'  'fdr' (default) | 'bonferroni' | 'none'
%   'fdr_scope'   'contrast' (default) -- correct across parcels within each
%                 contrast/term | 'all' -- correct across all parcels x contrasts
%   'transform'   'auto' (default) | 'fisherz' | 'none'  (contrast cells)
%   'map_value'   't' (default) | 'beta' -- what to paint into the maps
%
% Misc
%   'use_parallel'  logical (default false)
%   'verbose'       logical (default true)
%
% Output
% ------
%   results struct:
%     .per_parcel_rsms  [nParcels x 1] array of rsm objects
%     .contrast_table   long table (contrast path): Parcel | Contrast | t | p | FDR_P | sig
%     .lme_table        long table (lme path): Parcel | Term | beta | se | t | p | FDR_P | sig
%     .maps             struct: maps.<name> = statistic_image per contrast/term
%     .atlas            the atlas used
%     .history          cellstr

% =========================================================================
% Parse
% =========================================================================
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('atlas',        [],      @(x) isa(x, 'atlas'));
p.addParameter('group_by',     '',      @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('subject_var',  'subject_id', @(x) ischar(x) || isstring(x));
p.addParameter('metric',       'correlation', @(x) ischar(x) || isstring(x));
p.addParameter('contrasts',    {},      @iscell);
p.addParameter('lme',          '',      @(x) ischar(x) || isstring(x));
p.addParameter('predictors',   {},      @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('interactions', {},      @iscell);
p.addParameter('three_way',    {},      @iscell);
p.addParameter('tail',         'both',  @(x) ischar(x) || isstring(x));
p.addParameter('correction',   'fdr',   @(x) ischar(x) || isstring(x));
p.addParameter('fdr_scope',    'contrast', @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'contrast','all'}));
p.addParameter('transform',    'auto',  @(x) ischar(x) || isstring(x));
p.addParameter('map_value',    't',     @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'t','beta'}));
p.addParameter('use_parallel', false,   @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose',      true,    @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;
fwd = opt; % keep for compute_rsm/rsa_lme forwarding

if isempty(opt.atlas)
    error('rsa_parcelwise:noAtlas', 'Pass an atlas object via ''atlas''.');
end
do_contrasts = ~isempty(opt.contrasts);
do_lme       = ~isempty(char(opt.lme));
if ~do_contrasts && ~do_lme
    error('rsa_parcelwise:noAnalysis', 'Pass ''contrasts'' and/or ''lme''.');
end

atlas_obj = opt.atlas;
labels    = atlas_obj.labels;
n_parcels = num_regions(atlas_obj);

results = struct();
results.atlas   = atlas_obj;
results.maps    = struct();
results.history = {sprintf('%s: rsa_parcelwise on %d parcels', datestr(now), n_parcels)};

% =========================================================================
% Build per-parcel RSMs (reuse compute_rsm parcellation)
% =========================================================================
if opt.verbose, fprintf('rsa_parcelwise: building %d per-parcel RSMs...\n', n_parcels); end
extra = struct2nv(p.Unmatched);
R_pp = compute_rsm(dat, 'parcellation', atlas_obj, ...
    'group_by',    opt.group_by, ...
    'subject_var', opt.subject_var, ...
    'metric',      opt.metric, ...
    'verbose',     false, extra{:});
results.per_parcel_rsms = R_pp;

% =========================================================================
% CONTRAST PATH
% =========================================================================
if do_contrasts
    if opt.verbose, fprintf('rsa_parcelwise: running %d contrasts per parcel...\n', size(opt.contrasts,1)); end
    results.contrast_table = run_contrast_path(R_pp, labels, opt);
    results = build_maps_from_table(results, results.contrast_table, ...
        'Contrast', atlas_obj, opt);
end

% =========================================================================
% LME PATH
% =========================================================================
if do_lme
    if isempty(opt.predictors)
        error('rsa_parcelwise:noPredictors', 'LME path requires ''predictors''.');
    end
    if opt.verbose, fprintf('rsa_parcelwise: running per-parcel LME...\n'); end
    results.lme_table = run_lme_path(dat, atlas_obj, labels, opt);
    results = build_maps_from_table(results, results.lme_table, ...
        'Term', atlas_obj, opt);
end

if opt.verbose, fprintf('rsa_parcelwise: done. %d map(s) emitted.\n', numel(fieldnames(results.maps))); end

end


% =========================================================================
function T = run_contrast_path(R_pp, labels, opt)
% Run ttest_contrasts on each parcel; accumulate long table; FDR across parcels.

n_parcels = numel(R_pp);
spec = opt.contrasts;
n_con = size(spec, 1);
con_names = spec(:, 1);

rows = cell(n_parcels, 1);
for i = 1:n_parcels
    R = R_pp(i);
    if isempty(R) || size(R, 3) < 2 || isempty(fieldnames(R.groupings))
        continue
    end
    try
        Tc = R.ttest_contrasts(spec, 'tail', opt.tail, 'correction', 'none', ...
                               'transform', opt.transform);
    catch
        continue
    end
    parcel_col = repmat(labels(i), height(Tc), 1);
    rows{i} = table(parcel_col, string(Tc.Contrast), Tc.Mean_Diff, Tc.t, Tc.P, ...
        'VariableNames', {'Parcel', 'Contrast', 'effect', 't', 'p'});
end
T = vertcat(rows{:});
if isempty(T), return; end

% FDR across parcels (per contrast or all)
T = apply_parcel_correction(T, 'Contrast', con_names, opt);
end


% =========================================================================
function T = run_lme_path(dat, atlas_obj, labels, opt)
% Per-parcel apply_mask + rsa_lme; accumulate fixed-effect coefficients.

n_parcels = num_regions(atlas_obj);
rows = cell(n_parcels, 1);

for i = 1:n_parcels
    try
        parcel_mask = select_atlas_subset(atlas_obj, labels(i));
        dat_i = apply_mask(dat, parcel_mask);
        if isempty(dat_i.dat) || size(dat_i.dat, 1) < 3, continue; end
        mdl_i = dat_i.rsa_lme(char(opt.lme), ...
            'predictors',   cellstr(opt.predictors), ...
            'interactions', opt.interactions, ...
            'three_way',    opt.three_way, ...
            'subject_var',  opt.subject_var, ...
            'metric',       opt.metric, ...
            'verbose',      false);
        ce = mdl_i.Coefficients;
        % Drop intercept row
        keep = ~strcmp(ce.Name, '(Intercept)');
        terms = ce.Name(keep);
        parcel_col = repmat(labels(i), numel(terms), 1);
        rows{i} = table(parcel_col, string(terms), ce.Estimate(keep), ce.SE(keep), ...
            ce.tStat(keep), ce.pValue(keep), ...
            'VariableNames', {'Parcel', 'Term', 'beta', 'se', 't', 'p'});
    catch ME
        if opt.verbose
            warning('rsa_parcelwise:lmeParcelFailed', 'Parcel %s LME failed: %s', labels{i}, ME.message);
        end
    end
end
T = vertcat(rows{:});
if isempty(T), return; end

term_names = unique(T.Term, 'stable');
T = apply_parcel_correction(T, 'Term', cellstr(term_names), opt);
end


% =========================================================================
function T = apply_parcel_correction(T, group_col, group_names, opt)
% Add FDR_P + sig columns, correcting across parcels within each
% contrast/term (fdr_scope='contrast') or across everything ('all').

n = height(T);
fdr_p = nan(n, 1);
sig   = false(n, 1);
correction = lower(char(opt.correction));

if strcmpi(opt.fdr_scope, 'all')
    groups_to_do = {':'};
else
    groups_to_do = group_names(:)';
end

for g = 1:numel(groups_to_do)
    if strcmp(groups_to_do{g}, ':')
        idx = (1:n)';
    else
        idx = find(string(T.(group_col)) == string(groups_to_do{g}));
    end
    if isempty(idx), continue; end
    ps = T.p(idx);
    [adj, s] = correct_pvals(ps, correction);
    fdr_p(idx) = adj;
    sig(idx)   = s;
end

T.FDR_P = fdr_p;
T.sig   = sig;
end


% =========================================================================
function [adj, sig] = correct_pvals(ps, correction)
n = numel(ps);
switch correction
    case 'fdr'
        ranks = zeros(n,1); [~, order] = sort(ps); ranks(order) = 1:n;
        adj = min(ps(:) .* n ./ ranks, 1);
        if exist('FDR', 'file') == 2
            try
                thr = FDR(ps, 0.05);
                if isempty(thr) || isnan(thr), thr = 0; end
                sig = ps <= thr;
            catch
                sig = adj < 0.05;
            end
        else
            sig = adj < 0.05;
        end
    case 'bonferroni'
        adj = min(ps(:) * n, 1);
        sig = adj < 0.05;
    case 'none'
        adj = ps(:);
        sig = ps < 0.05;
    otherwise
        adj = ps(:); sig = ps < 0.05;
end
end


% =========================================================================
function results = build_maps_from_table(results, T, group_col, atlas_obj, opt)
% For each unique contrast/term, paint a statistic_image from the per-parcel
% values and store in results.maps.<name>.

if isempty(T), return; end
names = unique(string(T.(group_col)), 'stable');
val_col = 't';
if strcmpi(opt.map_value, 'beta') && ismember('beta', T.Properties.VariableNames)
    val_col = 'beta';
end

for i = 1:numel(names)
    rows = string(T.(group_col)) == names(i);
    parcels = cellstr(T.Parcel(rows));
    vals    = T.(val_col)(rows);
    pvals   = T.p(rows);

    map = assign_vals_to_atlas(atlas_obj, parcels, vals, ...
        'p_vals', pvals, 'output_type', 'statistic_image', ...
        'dat_descrip', sprintf('rsa_parcelwise %s = %s', group_col, names(i)));

    fld = matlab.lang.makeValidName(char(names(i)));
    results.maps.(fld) = map;
end
end


% =========================================================================
function nv = struct2nv(s)
% Convert an unmatched-args struct to a name-value cell for forwarding
fn = fieldnames(s);
nv = cell(1, 2*numel(fn));
for i = 1:numel(fn), nv{2*i-1} = fn{i}; nv{2*i} = s.(fn{i}); end
end
