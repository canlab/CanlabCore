function status = hrf_score_one_prefix(output_prefix, varargin)
%HRF_SCORE_ONE_PREFIX Score one (task, model) prefix's whole-brain HRF maps.
%
% status = hrf_score_one_prefix(output_prefix, ...)
%
% Single source of truth for applying signature/image sets to one set of
% whole-brain HRF maps. Writes *_<object>_map_scores.csv files in place at
% OUTPUT_PREFIX. Three callers share this helper:
%
%   * hrf_write_slurm_study_script (SLURM worker) - in-memory stats path,
%     called immediately after hrf_fit_wholebrain_stats with the in-memory
%     struct so no NIfTI re-read is needed.
%   * hrf_audit_slurm_outputs (RepairMissing mode) - disk path, called when
%     the audit detects core_complete && ~score_complete rows.
%   * hrf_score_wholebrain_input_table (post-hoc backfill) - disk path,
%     iterated per row of the second-level input table.
%
% Required input
%   output_prefix : char/string. Full path prefix used by
%                   hrf_fit_wholebrain_stats. The whole-brain image files
%                   are expected at <output_prefix>_beta.nii etc.
%
% Common name-value parameters
%   'ModelName'        - 'fir' | 'sfir' | 'canonical' | 'spline' | ''.
%                        Used to build per-model filenames when NumModels>1.
%   'NumModels'        - 1 if only one whole-brain model was fit at this
%                        prefix; >1 if multiple share the prefix. Controls
%                        filename composition: single-model files are
%                        <prefix>_<object>_map_scores.csv; multi-model files
%                        are <prefix>_<model>_<object>_map_scores.csv.
%   'ScoreObjects'     - cellstr/string. Default {'beta'}. Subset of
%                        {'beta','t'}.
%   'SignatureSets'    - cellstr/string. Passed to
%                        hrf_apply_maps_to_wholebrain.
%   'ImageSets'        - cellstr/string/image_vector. Passed to
%                        hrf_apply_maps_to_wholebrain.
%   'SimilarityMetric' - default 'dotproduct'.
%   'PropagateSE'      - default true. Beta-only and metric-restricted, same
%                        semantics as hrf_apply_maps_to_wholebrain.
%
% Metadata resolution (in order of preference)
%   'MetadataTable'    - pre-resolved metadata table; skip resolution.
%   'MetadataFile'     - explicit metadata CSV; loaded if exists.
%   'ResultMatFile'    - result.mat path; fallback if no MetadataFile.
%   Otherwise: <prefix>_metadata.csv, then <prefix>_<model>_metadata.csv
%   siblings, then return empty (and skip if RequireMetadata is true).
%
% Fast path (skip NIfTI re-read)
%   'StatsInput'       - struct returned by hrf_fit_wholebrain_stats (fields
%                        .b, .t, .metadata_table). If supplied, no disk
%                        reads are needed for the score/SE objects.
%
% Other controls
%   'Overwrite'        - default false. Force regeneration of existing CSVs.
%   'OverwriteStale'   - default true. Regenerate CSVs that exist but fail
%                        metadata/signature-coverage validation.
%   'RequireMetadata'  - default true. Skip writing if no metadata could be
%                        resolved.
%   'NoVerbose'        - default true. Suppress fmri_data load messages.
%   'WarningContext'   - propagated to hrf_apply_maps_to_wholebrain.
%
% Return
%   status struct with fields
%     .output_prefix       - char
%     .model_name          - char
%     .metadata_file       - char (or '')
%     .wrote_files         - cellstr of CSV paths produced this call
%     .skipped_existing    - cellstr of CSV paths left alone (already valid)
%     .skipped_stale       - cellstr of CSV paths left alone (stale, with
%                            OverwriteStale=false)
%     .errors              - struct array with .object and .message fields
%     .core_inputs_present - logical; true iff *_beta.nii exists (or .b is
%                            provided via StatsInput)
%
% See also: hrf_apply_maps_to_wholebrain, hrf_audit_slurm_outputs,
%           hrf_score_wholebrain_input_table, hrf_write_slurm_study_script.

p = inputParser;
p.addRequired('output_prefix', @(x) ischar(x) || isstring(x));
p.addParameter('ModelName', '', @(x) ischar(x) || isstring(x));
p.addParameter('NumModels', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('ScoreObjects', {'beta'}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('SignatureSets', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ImageSets', {}, @(x) ischar(x) || iscell(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('AtlasObj', [], @(x) isempty(x) || isa(x, 'atlas') || isa(x, 'image_vector'));
p.addParameter('AtlasName', '', @(x) ischar(x) || isstring(x));
p.addParameter('Regions', {}, @(x) iscell(x) || isstring(x) || ischar(x));
% Multi-atlas: score several atlases in one pass (columns namespaced by name,
% so give each a DISTINCT AtlasNames entry to avoid collisions). AtlasRegions
% is a cell-of-region-lists, one per atlas ({} = all regions for that atlas).
% The singular AtlasObj/AtlasName/Regions above still work and are scored too.
p.addParameter('AtlasObjs', {}, @(x) isempty(x) || iscell(x) || isa(x, 'atlas') || isa(x, 'image_vector'));
p.addParameter('AtlasNames', {}, @(x) isempty(x) || iscell(x) || ischar(x) || isstring(x));
p.addParameter('AtlasRegions', {}, @(x) isempty(x) || iscell(x));
p.addParameter('Normalize', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.addParameter('PropagateSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('SEScoreSuffix', '_se', @(x) ischar(x) || isstring(x));
p.addParameter('MetadataTable', table(), @(x) isempty(x) || istable(x));
p.addParameter('MetadataFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('ResultMatFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('StatsInput', [], @(x) isempty(x) || isstruct(x));
p.addParameter('Overwrite', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('OverwriteStale', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Append', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('RequireMetadata', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('NoVerbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('WarningContext', '', @(x) ischar(x) || isstring(x));
p.parse(output_prefix, varargin{:});
opts = p.Results;

prefix = char(opts.output_prefix);
model_name = char(opts.ModelName);
n_models = double(opts.NumModels);
objects = local_resolve_score_objects(opts.ScoreObjects);

status = struct( ...
    'output_prefix', prefix, ...
    'model_name', model_name, ...
    'metadata_file', '', ...
    'wrote_files', {{}}, ...
    'skipped_existing', {{}}, ...
    'skipped_stale', {{}}, ...
    'errors', struct('object', {}, 'message', {}), ...
    'core_inputs_present', false);

% Resolve metadata first; many downstream checks need it.
[metadata_table, metadata_file] = local_resolve_metadata( ...
    prefix, model_name, opts.MetadataTable, opts.MetadataFile, ...
    opts.ResultMatFile, opts.StatsInput);
status.metadata_file = metadata_file;

if isempty(metadata_table) && logical(opts.RequireMetadata)
    status.errors(end + 1) = struct('object', '', ...
        'message', sprintf('Cannot score %s: no metadata table resolvable from prefix, result.mat, or sibling CSVs.', prefix));
    return
end

% Check core inputs (cheaply, before iterating objects).
status.core_inputs_present = local_core_inputs_present(prefix, model_name, n_models, opts.StatsInput);
if ~status.core_inputs_present
    expected_beta = [local_image_prefix(prefix, model_name, n_models) '_beta.nii'];
    status.errors(end + 1) = struct('object', '', ...
        'message', sprintf('Missing core whole-brain inputs for prefix %s (expected %s or StatsInput.b).', prefix, expected_beta));
    return
end

for j = 1:numel(objects)
    object_name = objects{j};
    score_file = local_score_file_path(prefix, model_name, n_models, object_name);
    require_uncertainty = local_requires_uncertainty(prefix, object_name, opts);

    if exist(score_file, 'file') == 2 && ~logical(opts.Overwrite)
        if local_existing_score_is_valid(score_file, metadata_table, require_uncertainty, opts)
            status.skipped_existing{end + 1} = score_file; %#ok<AGROW>
            continue
        end
        % Existing file is "stale" relative to current request -- either
        % missing requested signature/image sets, or out-of-shape metadata.
        % If Append mode is on and metadata aligns, compute ONLY the missing
        % sets and merge into the existing CSV (preserves prior columns).
        % Otherwise fall back to OverwriteStale (full regenerate) or skip.
        if logical(opts.Append)
            [appended, append_err] = local_try_append( ...
                score_file, prefix, object_name, model_name, metadata_table, opts);
            if appended
                status.wrote_files{end + 1} = score_file; %#ok<AGROW>
                continue
            elseif ~isempty(append_err)
                status.errors(end + 1) = struct('object', object_name, ...
                    'message', sprintf('append failed for %s: %s (falling back to overwrite if OverwriteStale)', ...
                        score_file, append_err)); %#ok<AGROW>
                % fall through; OverwriteStale (if true) regenerates from scratch
            end
            % append_err empty AND not appended means metadata didn't align;
            % fall through to the OverwriteStale logic below.
        end
        if ~logical(opts.OverwriteStale)
            status.skipped_stale{end + 1} = score_file; %#ok<AGROW>
            continue
        end
        % else fall through and regenerate from scratch
    end

    try
        [score_obj, se_obj] = local_resolve_score_objects_for_write( ...
            prefix, object_name, metadata_table, opts);
        local_validate_score_object_metadata(score_obj, metadata_table, object_name, prefix);

        scores = hrf_apply_maps_to_wholebrain(score_obj, ...
            'Object', object_name, ...
            'SignatureSets', opts.SignatureSets, ...
            'ImageSets', opts.ImageSets, ...
            'SimilarityMetric', opts.SimilarityMetric, ...
            'SEInput', se_obj, ...
            'PropagateSE', opts.PropagateSE, ...
            'SEScoreSuffix', opts.SEScoreSuffix, ...
            'MetadataTable', metadata_table, ...
            'OutputCsv', '', ...
            'WarningContext', local_warning_context(opts.WarningContext, model_name, object_name, prefix));

        % Optionally append atlas region means -- one block per atlas spec.
        atlas_specs = local_atlas_specs(opts);
        for k = 1:numel(atlas_specs)
            atlas_cols = local_compute_atlas_scores(score_obj, atlas_specs(k).obj, ...
                atlas_specs(k).name, ...
                local_effective_regions(atlas_specs(k).obj, atlas_specs(k).regions), ...
                opts.Normalize);
            scores = local_horzappend_scores(scores, atlas_cols);
        end

        writetable(scores, score_file);
        status.wrote_files{end + 1} = score_file; %#ok<AGROW>
    catch err
        status.errors(end + 1) = struct('object', object_name, 'message', err.message); %#ok<AGROW>
    end
end
end


% =========================================================================
% Append mode: compute only missing sets, merge into existing CSV
% =========================================================================
function [appended, err_msg] = local_try_append(score_file, prefix, object_name, model_name, metadata_table, opts)
% Returns appended=true and err_msg='' on success.
% Returns appended=false and err_msg='' if metadata mismatch -> caller falls back.
% Returns appended=false and err_msg=<message> on actual error -> caller decides.
appended = false;
err_msg = '';

try
    existing = readtable(score_file, 'TextType', 'string');
catch err
    err_msg = sprintf('could not read existing CSV: %s', err.message);
    return
end

if ~local_existing_metadata_aligns(existing, metadata_table)
    % Silent fall-through; the existing file is incompatible with current
    % metadata (different rows, conditions, or lags). Full regenerate is
    % the right move.
    return
end

[missing_sigs, missing_imgs] = local_missing_sets(existing, opts);
atlas_specs = local_atlas_specs(opts);
missing_regions = local_missing_atlas_regions(existing, opts);
if isempty(missing_sigs) && isempty(missing_imgs) && isempty(missing_regions)
    % Everything requested is already present; nothing to do, but the
    % existing file is fine.
    appended = true;  % no-op write avoided
    return
end

try
    [score_obj, se_obj] = local_resolve_score_objects_for_write( ...
        prefix, object_name, metadata_table, opts);
    local_validate_score_object_metadata(score_obj, metadata_table, object_name, prefix);

    new_scores = hrf_apply_maps_to_wholebrain(score_obj, ...
        'Object', object_name, ...
        'SignatureSets', missing_sigs, ...
        'ImageSets', missing_imgs, ...
        'SimilarityMetric', opts.SimilarityMetric, ...
        'SEInput', se_obj, ...
        'PropagateSE', opts.PropagateSE, ...
        'SEScoreSuffix', opts.SEScoreSuffix, ...
        'MetadataTable', metadata_table, ...
        'OutputCsv', '', ...
        'WarningContext', local_warning_context(opts.WarningContext, model_name, object_name, prefix));

    % Atlas region means -- per atlas spec, append ONLY that atlas's missing
    % regions, so adding a new region (or a whole new atlas) to an existing
    % CSV doesn't trigger a full re-extract.
    for k = 1:numel(atlas_specs)
        mr = local_missing_atlas_regions_spec(existing, atlas_specs(k), opts.Normalize);
        if ~isempty(mr)
            atlas_cols = local_compute_atlas_scores(score_obj, atlas_specs(k).obj, ...
                atlas_specs(k).name, mr, opts.Normalize);
            new_scores = local_horzappend_scores(new_scores, atlas_cols);
        end
    end

    merged = local_merge_score_tables(existing, new_scores);
    writetable(merged, score_file);
    appended = true;
catch err
    err_msg = err.message;
end
end


function tf = local_existing_metadata_aligns(existing, meta)
% Row count + condition + lag must match for an in-place append to make sense.
tf = false;
if isempty(meta) || height(existing) ~= height(meta)
    return
end
has_cond_e = any(strcmp('condition', existing.Properties.VariableNames));
has_cond_m = any(strcmp('condition', meta.Properties.VariableNames));
if has_cond_e && has_cond_m && ~isequal(string(existing.condition), string(meta.condition))
    return
end
has_lag_e = any(strcmp('lag_index', existing.Properties.VariableNames));
has_lag_m = any(strcmp('lag_index', meta.Properties.VariableNames));
if has_lag_e && has_lag_m && ~isequaln(local_to_numeric(existing.lag_index), local_to_numeric(meta.lag_index))
    return
end
tf = true;
end


function [missing_sigs, missing_imgs] = local_missing_sets(existing, opts)
% Semantic match: each requested SET name produces columns starting with
% sig_<setname>_ or map_<setname>_. A set is "present" if at least one
% such column exists.
missing_sigs = {};
sigsets = local_to_cell(opts.SignatureSets);
for i = 1:numel(sigsets)
    if ~local_has_numeric_prefix(existing, local_var_prefix({'sig', sigsets{i}}))
        missing_sigs{end + 1} = sigsets{i}; %#ok<AGROW>
    end
end
missing_imgs = {};
image_sets = local_to_cell(opts.ImageSets);
for i = 1:numel(image_sets)
    if isa(image_sets{i}, 'image_vector')
        set_name = 'imageset';
    else
        set_name = char(image_sets{i});
    end
    if ~local_has_numeric_prefix(existing, local_var_prefix({'map', set_name}))
        missing_imgs{end + 1} = image_sets{i}; %#ok<AGROW>
    end
end
end


function merged = local_merge_score_tables(existing, new_scores)
% Add any new score columns from new_scores into existing. Metadata
% columns come from existing (they should match new_scores's metadata
% since the alignment check passed). Existing score columns are
% preserved -- if a column name collides we keep the existing one
% (assume it's already valid).
merged = existing;
metadata_cols = {'volume_index', 'condition', 'condition_index', 'lag_index', ...
    'lag_seconds', 'image_label', 'N', 'dfe', 'TR', 'mode', 'subject', 'run_label'};
new_vars = new_scores.Properties.VariableNames;
existing_vars = existing.Properties.VariableNames;
for i = 1:numel(new_vars)
    col = new_vars{i};
    if ismember(col, metadata_cols), continue; end
    if any(strcmp(col, existing_vars)), continue; end
    merged.(col) = new_scores.(col);
end
end


% =========================================================================
% Atlas region-mean scoring (post-hoc extraction from whole-brain maps)
% =========================================================================
function atlas_cols = local_compute_atlas_scores(score_obj, atlas_obj, atlas_name, regions, normalize)
% Returns a table with one column per requested atlas region:
%   atlas_<atlasname>_<region>_<suffix>
% with row count == n_volumes in score_obj.
% Uses extract_roi_averages which handles space resampling between atlas
% and score_obj automatically. If `regions` is non-empty, subsets the
% atlas via select_atlas_subset(..., 'exact', 'deterministic') first.
% The 'deterministic' flag forces probabilistic atlases (e.g., canlab2024)
% to a winner-takes-all parcellation so each voxel belongs to at most one
% region -- avoids the probabilistic overlap that would otherwise blur
% region boundaries.
%
% normalize is one of:
%   'mean'  - region mean (extract_roi_averages default). Suffix '_mean'.
%             Already L1-normalized by region size (= sum / nVoxels).
%             (Default.) Standard ROI output: plain mean BOLD signal.
%   'l1'    - region mean, then each region's time course divided by its
%             own L1 norm sum(|values|). Makes regions of differing baseline
%             signal magnitudes directly comparable in HRF *shape*. Suffix
%             '_meanL1'.
%   'none'  - region sum (extract_roi_averages 'nonorm'). Suffix '_sum'.
%             Larger regions get larger values; only useful when downstream
%             code wants raw integrated signal.
if nargin < 5 || isempty(normalize), normalize = 'mean'; end
normalize = lower(strtrim(char(normalize)));

if isempty(regions)
    regions = local_effective_regions(atlas_obj, {});
end
regions = cellstr(string(regions));

% IMPORTANT: extract_roi_averages exhibits different per-volume behavior
% when given a multi-region atlas (518 regions in one call) vs a
% single-region deterministic subset. The multi-region path can produce
% jagged per-volume averages even when each region's underlying voxels
% are smooth canonical curves -- confirmed in a side-by-side diagnostic
% where extract_roi_averages on a 1-region subset gave a smooth peak
% while the CSV column (written from a 518-region subset call) was
% noise. To stay on the correct code path, subset the atlas to ONE
% region at a time and call extract_roi_averages per region. Slower
% (~Nx for N regions) but correct.

extract_args = {};
switch normalize
    case 'none'
        extract_args = {'nonorm'};
        suffix = 'sum';
    case {'mean', 'l1'}
        suffix = normalize_suffix(normalize);
    otherwise
        error('hrf_score_one_prefix:UnknownNormalize', ...
            'Unknown Normalize: %s. Use mean, l1, or none.', normalize);
end

atlas_token = matlab.lang.makeValidName(char(atlas_name));
n_vol = size(score_obj.dat, 2);
atlas_cols = table();

% Convert score_obj to fmri_data for extract_roi_averages. The scoring
% pipeline loads NIfTIs as statistic_image(fmri_data(...)); on that
% statistic_image-typed input extract_roi_averages takes a different
% code path that applies sig() filters and produces non-linear per-volume
% behavior even for smooth canonical data. The user's manual diagnostic
% on plain fmri_data gave smooth output; matching that here.
extract_input = local_to_fmri_data_for_extract(score_obj);

for r = 1:numel(regions)
    region_label = regions{r};
    if isa(atlas_obj, 'atlas')
        % Always 'deterministic' (winner-take-all). Probabilistic atlases are
        % otherwise far too liberal -- a canlab2024 region spans ~2x the voxels
        % of its hard parcellation, blurring boundaries past what is canonically
        % anatomical. On a label-only atlas (empty probability_maps) the flag is
        % a verified harmless no-op (it falls back to the integer .dat labels).
        try
            single_sub = select_atlas_subset(atlas_obj, {region_label}, 'exact', 'deterministic');
        catch err
            warning('hrf_score_one_prefix:AtlasRegionSubset', ...
                'Skipping region ''%s'': %s', region_label, err.message);
            continue
        end
    else
        single_sub = atlas_obj;
    end

    try
        cl = extract_roi_averages(extract_input, single_sub, extract_args{:});
    catch err
        warning('hrf_score_one_prefix:AtlasExtract', ...
            'Skipping region ''%s'': %s', region_label, err.message);
        continue
    end
    if isempty(cl)
        continue
    end

    region_token = matlab.lang.makeValidName(char(region_label));
    col_name = local_make_atlas_column_name(atlas_token, region_token, suffix);
    vals = double(cl(1).dat(:));
    if numel(vals) ~= n_vol
        padded = NaN(n_vol, 1);
        m = min(numel(vals), n_vol);
        padded(1:m) = vals(1:m);
        vals = padded;
    end
    if strcmp(normalize, 'l1')
        l1 = sum(abs(vals), 'omitnan');
        if l1 > 0
            vals = vals / l1;
        end
    end
    atlas_cols.(col_name) = vals;
end
end


function fd = local_to_fmri_data_for_extract(obj)
% Return a plain fmri_data so extract_roi_averages takes the fmri_data code
% path (a statistic_image triggers a per-volume mask path that gives jagged
% atlas output). The conversion MUST yield a valid object whose mask/volInfo
% let extract_roi_averages resample the atlas to the data space.
%
% Use the fmri_data(obj) constructor: it builds a valid object directly from
% the IN-MEMORY data (no disk reload), so it works even for a StatsInput score
% object that has no fullpath -- and it preserves the raw .dat (the .sig
% threshold mask of a statistic_image is NOT applied), matching the historical
% reload-from-disk behavior. Verified against a real 2mm-atlas vs 2.683mm-beta
% mismatch: the previous hand field-copy raised "Arrays have incompatible
% sizes" on the resample (every region skipped, no atlas columns written),
% while fmri_data(obj) succeeds. Disk reload / field-copy are kept only as
% last-resort fallbacks.
if isa(obj, 'fmri_data') && ~isa(obj, 'statistic_image')
    fd = obj;
    return
end

% Primary: convert in memory (handles statistic_image / image_vector).
try
    fd = fmri_data(obj);
    return
catch
    % fall through to disk reload / field copy
end

% Fallback 1: reload from disk if the object is file-backed.
if isprop(obj, 'fullpath') && ~isempty(obj.fullpath)
    src = deblank(obj.fullpath(1, :));
    if exist(src, 'file') == 2
        try
            fd = fmri_data(src, 'noverbose');
            return
        catch
        end
    end
end

% Fallback 2 (last resort): copy shared fields. Known to break the atlas
% resample, but better than nothing if both conversions above failed.
fd = fmri_data();
shared = {'dat', 'volInfo', 'removed_voxels', 'removed_images', ...
    'space_defining_image_name', 'fullpath', 'files_exist', 'history', ...
    'image_names', 'source_notes', 'mask', 'mask_descrip', ...
    'images_per_session'};
for i = 1:numel(shared)
    f = shared{i};
    if isprop(fd, f) && isprop(obj, f)
        try
            fd.(f) = obj.(f);
        catch
        end
    end
end
end


function s = normalize_suffix(mode)
switch lower(char(mode))
    case 'mean', s = 'mean';
    case 'l1',   s = 'meanL1';
    case 'none', s = 'sum';
    otherwise,   s = lower(char(mode));
end
end


function regions = local_effective_regions(atlas_obj, requested_regions)
% Returns the cellstr of region labels to score:
%   - if requested_regions is non-empty -> that list (cellstr)
%   - else atlas_obj.labels (full atlas)
%   - else generic region_NNN fallback
if ~isempty(requested_regions)
    regions = cellstr(string(requested_regions));
    regions = regions(:)';
    return
end
if isprop(atlas_obj, 'labels') && ~isempty(atlas_obj.labels)
    regions = cellstr(string(atlas_obj.labels));
    regions = regions(:)';
    return
end
% Last-ditch fallback: ask the atlas for a region count if possible.
n_guess = 0;
if isprop(atlas_obj, 'dat') && ~isempty(atlas_obj.dat)
    n_guess = max(max(double(atlas_obj.dat(:))), 0);
end
if n_guess <= 0
    n_guess = 1;
end
regions = arrayfun(@(i) sprintf('region_%03d', i), 1:n_guess, 'UniformOutput', false);
end


function labels = local_atlas_labels(atlas_obj, n_regions)
labels = {};
if isprop(atlas_obj, 'labels') && ~isempty(atlas_obj.labels)
    labels = cellstr(string(atlas_obj.labels));
end
if numel(labels) < n_regions
    extra = arrayfun(@(i) sprintf('region_%03d', i), numel(labels) + 1 : n_regions, ...
        'UniformOutput', false);
    labels = [labels(:); extra(:)];
end
labels = labels(1:n_regions);
end


function name = local_resolve_atlas_name(atlas_obj, user_name)
% Priority: explicit user-supplied name -> atlas_obj.atlas_name property
% -> 'atlas' fallback.
name = char(user_name);
if ~isempty(name), return; end
if isprop(atlas_obj, 'atlas_name') && ~isempty(atlas_obj.atlas_name)
    name = char(atlas_obj.atlas_name);
    return
end
name = 'atlas';
end


function missing_regions = local_missing_atlas_regions(existing, opts)
% Aggregate missing-region labels across ALL atlas specs (used by the
% validity/no-op checks, which only need to know whether anything is
% missing). Per-atlas append uses local_missing_atlas_regions_spec.
missing_regions = {};
specs = local_atlas_specs(opts);
for k = 1:numel(specs)
    mr = local_missing_atlas_regions_spec(existing, specs(k), opts.Normalize);
    missing_regions = [missing_regions, mr]; %#ok<AGROW>
end
end

function mr = local_missing_atlas_regions_spec(existing, spec, normalize)
% Region labels for one atlas whose atlas_<name>_<region>_<suffix> column is
% absent from the existing CSV (suffix depends on Normalize). Different
% normalization modes produce different suffixes, so a CSV with _mean columns
% and a request for _meanL1 correctly routes the L1 versions through append.
mr = {};
atlas_token = matlab.lang.makeValidName(char(spec.name));
effective_regions = local_effective_regions(spec.obj, spec.regions);
suffix = normalize_suffix(normalize);
existing_vars = existing.Properties.VariableNames;
for r = 1:numel(effective_regions)
    region_token = matlab.lang.makeValidName(char(effective_regions{r}));
    col_name = local_make_atlas_column_name(atlas_token, region_token, suffix);
    if ~any(strcmp(col_name, existing_vars))
        mr{end + 1} = effective_regions{r}; %#ok<AGROW>
    end
end
end

function specs = local_atlas_specs(opts)
% Unify the singular (AtlasObj/AtlasName/Regions) and plural (AtlasObjs/
% AtlasNames/AtlasRegions) atlas inputs into a struct array of specs, each
% with fields .obj, .name, .regions. The singular spec (if any) comes first.
specs = struct('obj', {}, 'name', {}, 'regions', {});

if ~isempty(opts.AtlasObj)
    specs = local_append_atlas_spec(specs, opts.AtlasObj, ...
        local_resolve_atlas_name(opts.AtlasObj, opts.AtlasName), local_to_cell(opts.Regions));
end

objs = opts.AtlasObjs;
if ~isempty(objs)
    if ~iscell(objs), objs = {objs}; end
    names = opts.AtlasNames;
    if ischar(names) || isstring(names), names = cellstr(string(names)); end
    if ~iscell(names), names = {}; end
    regs = opts.AtlasRegions;
    if ~iscell(regs), regs = {}; end
    for k = 1:numel(objs)
        if isempty(objs{k}), continue; end
        if k <= numel(names) && ~isempty(names{k})
            nm = char(string(names{k}));
        else
            nm = local_resolve_atlas_name(objs{k}, '');
        end
        if k <= numel(regs), rg = local_to_cell(regs{k}); else, rg = {}; end
        specs = local_append_atlas_spec(specs, objs{k}, nm, rg);
    end
end

% Robustness: two atlases with the same resolved name would write colliding
% atlas_<name>_<region> columns (common when both inherit the same internal
% atlas_name). Auto-disambiguate later duplicates with a numeric suffix.
specs = local_disambiguate_spec_names(specs);
end

function specs = local_disambiguate_spec_names(specs)
seen = {};
for k = 1:numel(specs)
    base = char(specs(k).name);
    token = matlab.lang.makeValidName(base);
    if any(strcmp(token, seen))
        n = 2;
        % Counter PREFIX (not suffix): the 63-char column-name cap trims the
        % END of the atlas token, so a trailing suffix could be cut and collide
        % again; a prefix always survives.
        while any(strcmp(matlab.lang.makeValidName(sprintf('dup%d_%s', n, base)), seen))
            n = n + 1;
        end
        new_name = sprintf('dup%d_%s', n, base);
        warning('hrf_score_one_prefix:DuplicateAtlasName', ...
            ['Two atlases resolved to the same name ''%s''; renaming the later one to ''%s'' so ' ...
             'their region columns do not collide. Pass distinct AtlasNames to control this.'], ...
            base, new_name);
        specs(k).name = new_name;
        token = matlab.lang.makeValidName(new_name);
    end
    seen{end + 1} = token; %#ok<AGROW>
end
end

function specs = local_append_atlas_spec(specs, obj, name, regions)
n = numel(specs) + 1;
specs(n).obj = obj;
specs(n).name = name;
specs(n).regions = regions;
end


function col_name = local_make_atlas_column_name(atlas_token, region_token, suffix)
% Build 'atlas_<atlas_token>_<region_token>_<suffix>' but keep it under
% the historical MATLAB name length limit (63 chars) so the column can
% be assigned to a table variable. We hard-cap at 63 even on newer
% MATLAB releases where namelengthmax is 2048 -- this keeps column names
% identical across local dev and cluster (where namelengthmax may
% differ), so append/missing-region detection works across systems.
% The atlas_token is shortened if needed; the region_token is preserved
% at full length because it's the more informative half. If the region
% alone exceeds the budget, it's truncated as a last resort.
overhead = numel('atlas_') + 1 + 1 + numel(suffix);  % 'atlas_' + '_' + '_' + suffix
max_len = min(63, namelengthmax);
budget = max_len - overhead;

if numel(atlas_token) + numel(region_token) > budget
    atlas_budget = max(1, budget - numel(region_token));
    if numel(atlas_token) > atlas_budget
        atlas_token = atlas_token(1:atlas_budget);
    end
end

if numel(atlas_token) + numel(region_token) > budget
    region_budget = max(1, budget - numel(atlas_token));
    if numel(region_token) > region_budget
        region_token = region_token(1:region_budget);
    end
end

col_name = sprintf('atlas_%s_%s_%s', atlas_token, region_token, suffix);
end


function merged = local_horzappend_scores(scores, atlas_cols)
% Append atlas_cols as new columns to scores; row counts must match.
if isempty(atlas_cols) || width(atlas_cols) == 0
    merged = scores;
    return
end
if height(scores) ~= height(atlas_cols)
    error('hrf_score_one_prefix:AtlasRowMismatch', ...
        'Atlas scores have %d rows but main scores have %d rows.', ...
        height(atlas_cols), height(scores));
end
existing_vars = scores.Properties.VariableNames;
merged = scores;
new_vars = atlas_cols.Properties.VariableNames;
for i = 1:numel(new_vars)
    col = new_vars{i};
    if any(strcmp(col, existing_vars)), continue; end
    merged.(col) = atlas_cols.(col);
end
end

% =========================================================================
% Score-objects resolution (validate and normalize ScoreObjects argument)
% =========================================================================
function objects = local_resolve_score_objects(score_objects)
objects = lower(cellstr(string(score_objects)));
objects = cellfun(@strtrim, objects, 'UniformOutput', false);
if any(strcmp(objects, 'both')) || any(strcmp(objects, 'all'))
    objects = {'beta', 't'};
end
objects = unique(objects(~cellfun(@isempty, objects)), 'stable');
valid = {'beta', 'b', 't', 'tmap', 'tmaps'};
bad = setdiff(objects, valid);
if ~isempty(bad)
    error('hrf_score_one_prefix:UnknownScoreObjects', ...
        'Unknown ScoreObjects: %s. Use beta, t, or both.', strjoin(bad, ', '));
end
objects(strcmp(objects, 'b')) = {'beta'};
objects(ismember(objects, {'tmap', 'tmaps'})) = {'t'};
objects = unique(objects, 'stable');
end

% =========================================================================
% Filename composition
%   The helper's prefix arg is the task-level output_prefix. For multi-model
%   studies, the per-model NIfTI files live at <prefix>_<model>_*.nii, and
%   the per-model score CSV at <prefix>_<model>_<object>_map_scores.csv.
%   For single-model studies the model suffix is omitted from both.
% =========================================================================
function img_prefix = local_image_prefix(prefix, model_name, n_models)
if n_models <= 1 || isempty(model_name)
    img_prefix = prefix;
else
    img_prefix = sprintf('%s_%s', prefix, lower(char(model_name)));
end
end

function score_file = local_score_file_path(prefix, model_name, n_models, object_name)
if n_models <= 1 || isempty(model_name)
    score_file = sprintf('%s_%s_map_scores.csv', prefix, object_name);
else
    score_file = sprintf('%s_%s_%s_map_scores.csv', prefix, lower(model_name), object_name);
end
end

% =========================================================================
% Core-inputs presence check (cheap)
% =========================================================================
function tf = local_core_inputs_present(prefix, model_name, n_models, stats_input)
if ~isempty(stats_input) && isstruct(stats_input) && isfield(stats_input, 'b') && ...
        isobject(stats_input.b) && ~isempty(stats_input.b.dat)
    tf = true;
    return
end
img_prefix = local_image_prefix(prefix, model_name, n_models);
tf = exist([img_prefix '_beta.nii'], 'file') == 2;
end

% =========================================================================
% Metadata resolution
%   Order: MetadataTable input -> StatsInput.metadata_table ->
%          MetadataFile input -> <prefix>_metadata.csv ->
%          result.mat fallback -> sibling-model CSV fallback.
% =========================================================================
function [metadata_table, metadata_file] = local_resolve_metadata( ...
    prefix, model_name, metadata_in, metadata_file_in, result_mat_in, stats_input)

metadata_table = table();
metadata_file = '';

if ~isempty(metadata_in) && istable(metadata_in) && height(metadata_in) > 0
    metadata_table = metadata_in;
    if ~isempty(metadata_file_in)
        metadata_file = char(metadata_file_in);
    end
    return
end

if ~isempty(stats_input) && isstruct(stats_input) && ...
        isfield(stats_input, 'metadata_table') && istable(stats_input.metadata_table) && ...
        height(stats_input.metadata_table) > 0
    metadata_table = stats_input.metadata_table;
    return
end

if ~isempty(metadata_file_in)
    metadata_file = char(metadata_file_in);
    if exist(metadata_file, 'file') == 2
        try
            metadata_table = readtable(metadata_file, 'TextType', 'string');
            return
        catch
        end
    end
end

% Default: <prefix>_metadata.csv
default_metadata_file = [char(prefix) '_metadata.csv'];
if exist(default_metadata_file, 'file') == 2
    try
        metadata_table = readtable(default_metadata_file, 'TextType', 'string');
        metadata_file = default_metadata_file;
        return
    catch
    end
end

% Result MAT fallback
result_mat = char(result_mat_in);
if ~isempty(result_mat) && exist(result_mat, 'file') == 2
    metadata_table = local_metadata_from_result_mat(result_mat, model_name);
    if ~isempty(metadata_table)
        if isempty(metadata_file), metadata_file = default_metadata_file; end
        try
            writetable(metadata_table, metadata_file);
        catch
            metadata_file = '';
        end
        return
    end
end

% Sibling-model CSV fallback
sibling_metadata = local_metadata_from_sibling(prefix, model_name);
if ~isempty(sibling_metadata)
    metadata_table = sibling_metadata;
    if isempty(metadata_file), metadata_file = default_metadata_file; end
    try
        writetable(metadata_table, metadata_file);
    catch
        metadata_file = '';
    end
end
end

function metadata_table = local_metadata_from_result_mat(mat_file, model_name)
metadata_table = table();
try
    S = load(mat_file, 'results');
catch
    return
end
if ~isfield(S, 'results'), return; end
R = S.results;
if isempty(model_name), model_name = 'fir'; end
model_field = matlab.lang.makeValidName(lower(char(model_name)));

if isfield(R, 'wholebrain_metadata_by_model') && isfield(R.wholebrain_metadata_by_model, model_field)
    metadata_table = R.wholebrain_metadata_by_model.(model_field);
elseif isfield(R, 'wholebrain_by_model') && isfield(R.wholebrain_by_model, model_field) && ...
        isfield(R.wholebrain_by_model.(model_field), 'metadata_table')
    metadata_table = R.wholebrain_by_model.(model_field).metadata_table;
elseif isfield(R, 'wholebrain_metadata_table')
    metadata_table = R.wholebrain_metadata_table;
elseif isfield(R, 'wholebrain') && isfield(R.wholebrain, 'metadata_table')
    metadata_table = R.wholebrain.metadata_table;
end
end

function metadata_table = local_metadata_from_sibling(prefix, model_name)
metadata_table = table();
base_prefix = local_base_prefix(prefix, model_name);
candidate_files = { ...
    [base_prefix '_metadata.csv'], ...
    [base_prefix '_sfir_metadata.csv'], ...
    [base_prefix '_canonical_metadata.csv'], ...
    [base_prefix '_spline_metadata.csv']};
candidate_files = unique(candidate_files, 'stable');
this_file = [prefix '_metadata.csv'];

for i = 1:numel(candidate_files)
    fname = candidate_files{i};
    if strcmp(fname, this_file) || exist(fname, 'file') ~= 2
        continue
    end
    try
        T = readtable(fname, 'TextType', 'string');
    catch
        continue
    end
    if any(strcmp('condition', T.Properties.VariableNames)) && ...
            any(strcmp('lag_index', T.Properties.VariableNames))
        metadata_table = local_relabel_metadata(T, model_name);
        return
    end
end
end

function base_prefix = local_base_prefix(prefix, model_name)
base_prefix = char(prefix);
if isempty(model_name), return; end
model_suffix = ['_' lower(char(model_name))];
if endsWith(lower(base_prefix), model_suffix)
    base_prefix = base_prefix(1:end - numel(model_suffix));
end
end

function T = local_relabel_metadata(T, model_name)
model_name = lower(char(model_name));
T.volume_index = (1:height(T))';
if any(strcmp('mode', T.Properties.VariableNames))
    T.mode = repmat(string(upper(model_name)), height(T), 1);
end
if any(strcmp('image_label', T.Properties.VariableNames)) && ...
        any(strcmp('condition', T.Properties.VariableNames)) && ...
        any(strcmp('lag_index', T.Properties.VariableNames))
    T.image_label = local_image_labels(T, model_name);
end
end

function labels = local_image_labels(T, model_name)
conditions = cellstr(string(T.condition));
lag_index = local_to_numeric(T.lag_index);
if any(strcmp('lag_seconds', T.Properties.VariableNames))
    lag_seconds = local_to_numeric(T.lag_seconds);
else
    lag_seconds = lag_index;
end
labels = cell(height(T), 1);
for i = 1:height(T)
    if ismember(model_name, {'canonical', 'spline'})
        labels{i} = sprintf('%s_%s_lag%03d_%0.3gs', ...
            matlab.lang.makeValidName(model_name), ...
            matlab.lang.makeValidName(conditions{i}), ...
            lag_index(i), lag_seconds(i));
    else
        labels{i} = sprintf('%s_lag%03d_%0.3gs', ...
            matlab.lang.makeValidName(conditions{i}), lag_index(i), lag_seconds(i));
    end
end
end

% =========================================================================
% Existing-score validity
%   We regenerate when the existing CSV is missing requested signature/map
%   columns, when its row count or condition/lag layout disagrees with the
%   current metadata, or when uncertainty is requested but absent.
% =========================================================================
function tf = local_existing_score_is_valid(score_file, metadata_table, require_uncertainty, opts)
tf = true;
try
    S = readtable(score_file, 'TextType', 'string');
catch
    tf = false;
    return
end
if logical(require_uncertainty) && ~local_score_table_has_uncertainty(S)
    tf = false;
    return
end
if ~local_has_requested_score_sets(S, opts)
    tf = false;
    return
end
if isempty(metadata_table)
    return
end
if height(S) ~= height(metadata_table)
    tf = false;
    return
end
has_score_condition = any(strcmp('condition', S.Properties.VariableNames));
has_metadata_condition = any(strcmp('condition', metadata_table.Properties.VariableNames));
if has_metadata_condition && ~has_score_condition
    tf = false;
    return
elseif has_score_condition && has_metadata_condition && ...
        ~isequal(string(S.condition), string(metadata_table.condition))
    tf = false;
    return
end
has_score_lag = any(strcmp('lag_index', S.Properties.VariableNames));
has_metadata_lag = any(strcmp('lag_index', metadata_table.Properties.VariableNames));
if has_metadata_lag && ~has_score_lag
    tf = false;
    return
elseif has_score_lag && has_metadata_lag && ...
        ~isequaln(local_to_numeric(S.lag_index), local_to_numeric(metadata_table.lag_index))
    tf = false;
end
end

function tf = local_has_requested_score_sets(S, opts)
tf = true;
sigsets = local_to_cell(opts.SignatureSets);
for i = 1:numel(sigsets)
    if ~local_has_numeric_prefix(S, local_var_prefix({'sig', sigsets{i}}))
        tf = false;
        return
    end
end

image_sets = local_to_cell(opts.ImageSets);
for i = 1:numel(image_sets)
    if isa(image_sets{i}, 'image_vector')
        set_name = 'imageset';
    else
        set_name = char(image_sets{i});
    end
    if ~local_has_numeric_prefix(S, local_var_prefix({'map', set_name}))
        tf = false;
        return
    end
end

% Atlas region means: per-region check so adding a new region (or a new
% atlas) to an existing CSV doesn't require regenerating everything.
if ~isempty(local_missing_atlas_regions(S, opts))
    tf = false;
    return
end
end

function tf = local_has_numeric_prefix(S, prefix)
tf = false;
names = S.Properties.VariableNames;
for i = 1:numel(names)
    name = names{i};
    if startsWith(name, prefix) && isnumeric(S.(name)) && ~local_is_uncertainty_column(name)
        tf = true;
        return
    end
end
end

function prefix = local_var_prefix(parts)
prefix = matlab.lang.makeValidName(strjoin(cellfun(@char, parts, 'UniformOutput', false), '_'));
prefix = [prefix '_'];
end

function tf = local_score_table_has_uncertainty(S)
score_cols = local_score_columns_for_validation(S);
tf = true;
for i = 1:numel(score_cols)
    if ~any(strcmp([score_cols{i} '_se'], S.Properties.VariableNames))
        tf = false;
        return
    end
end
end

function cols = local_score_columns_for_validation(S)
metadata_cols = {'volume_index', 'condition', 'condition_index', 'lag_index', ...
    'lag_seconds', 'image_label', 'N', 'dfe', 'TR', 'mode'};
cols = {};
for i = 1:numel(S.Properties.VariableNames)
    name = S.Properties.VariableNames{i};
    if ismember(name, metadata_cols) || local_is_uncertainty_column(name)
        continue
    end
    if isnumeric(S.(name))
        cols{end + 1} = name; %#ok<AGROW>
    end
end
end

function tf = local_is_uncertainty_column(name)
tf = endsWith(char(name), '_se');
end

function tf = local_requires_uncertainty(prefix, object_name, opts)
img_prefix = local_image_prefix(prefix, char(opts.ModelName), double(opts.NumModels));
tf = strcmpi(char(object_name), 'beta') && logical(opts.PropagateSE) && ...
    local_is_linear_metric(opts.SimilarityMetric) && ...
    exist([img_prefix '_beta.nii'], 'file') == 2 && exist([img_prefix '_se.nii'], 'file') == 2;
end

function tf = local_is_linear_metric(metric)
metric = lower(strrep(strtrim(char(metric)), '_', ''));
tf = ismember(metric, {'dotproduct', 'dot'});
end

% =========================================================================
% Score / SE object loading
%   Fast path: take from StatsInput struct in memory (no disk read).
%   Slow path: load statistic_image(fmri_data(file)) from NIfTI.
% =========================================================================
function [score_obj, se_obj] = local_resolve_score_objects_for_write( ...
    prefix, object_name, metadata_table, opts)

if ~isempty(opts.StatsInput) && isstruct(opts.StatsInput)
    [score_obj, se_obj] = local_score_objects_from_stats(opts.StatsInput, object_name, opts);
    if ~isempty(score_obj)
        return
    end
end

% Disk path
model_name = char(opts.ModelName);
n_models = double(opts.NumModels);
score_obj = local_load_score_object_from_disk(prefix, object_name, metadata_table, ...
    logical(opts.NoVerbose), model_name, n_models);
se_obj = local_uncertainty_object_from_disk(prefix, object_name, metadata_table, ...
    logical(opts.NoVerbose), logical(opts.PropagateSE), model_name, n_models);
end

function [score_obj, se_obj] = local_score_objects_from_stats(stats_input, object_name, opts)
score_obj = [];
se_obj = [];
switch lower(char(object_name))
    case 'beta'
        if isfield(stats_input, 'b'), score_obj = stats_input.b; end
        if logical(opts.PropagateSE) && ~isempty(score_obj) && ...
                isprop(score_obj, 'ste') && ~isempty(score_obj.ste)
            se_obj = score_obj;
            se_obj.dat = score_obj.ste;
            se_obj.type = 'HRF beta standard error';
        end
    case 't'
        if isfield(stats_input, 't'), score_obj = stats_input.t; end
end
end

function obj = local_load_score_object_from_disk(prefix, object_name, metadata_table, noverbose, model_name, n_models)
image_file = local_score_image_file(prefix, object_name, model_name, n_models);
if exist(image_file, 'file') ~= 2
    error('hrf_score_one_prefix:MissingImage', ...
        'Missing %s image: %s', object_name, image_file);
end

load_args = {};
if noverbose, load_args = {'noverbose'}; end
obj = statistic_image(fmri_data(image_file, load_args{:}));
switch lower(char(object_name))
    case 'beta'
        obj.type = 'FIR HRF beta';
    case 't'
        obj.type = 'T';
end

obj = local_apply_metadata_to_obj(obj, metadata_table);
end

function obj = local_uncertainty_object_from_disk(prefix, object_name, metadata_table, noverbose, propagate_se, model_name, n_models)
obj = [];
if ~propagate_se || ~strcmpi(char(object_name), 'beta')
    return
end
image_file = [local_image_prefix(prefix, model_name, n_models) '_se.nii'];
if exist(image_file, 'file') ~= 2
    return
end

load_args = {};
if noverbose, load_args = {'noverbose'}; end
obj = statistic_image(fmri_data(image_file, load_args{:}));
obj.type = 'HRF beta standard error';
obj = local_apply_metadata_to_obj(obj, metadata_table);
end

function obj = local_apply_metadata_to_obj(obj, metadata_table)
if ~isempty(metadata_table) && height(metadata_table) == size(obj.dat, 2)
    if any(strcmp('image_label', metadata_table.Properties.VariableNames))
        obj.image_labels = cellstr(string(metadata_table.image_label));
    end
    if any(strcmp('N', metadata_table.Properties.VariableNames))
        obj.N = metadata_table.N(1);
    end
    if any(strcmp('dfe', metadata_table.Properties.VariableNames))
        obj.dfe = metadata_table.dfe(1);
    end
end
if isprop(obj, 'sig') && isempty(obj.sig)
    obj.sig = true(size(obj.dat));
end
end

function image_file = local_score_image_file(prefix, object_name, model_name, n_models)
img_prefix = local_image_prefix(prefix, model_name, n_models);
switch lower(char(object_name))
    case 'beta'
        image_file = [img_prefix '_beta.nii'];
    case 't'
        image_file = [img_prefix '_t.nii'];
    otherwise
        error('hrf_score_one_prefix:UnknownObject', ...
            'Unknown object: %s. Use beta or t.', char(object_name));
end
end

function local_validate_score_object_metadata(obj, metadata_table, object_name, prefix)
if isempty(metadata_table)
    return
end
n_images = size(obj.dat, 2);
if height(metadata_table) ~= n_images
    error('hrf_score_one_prefix:MetadataVolumeMismatch', ...
        'Cannot score %s: metadata has %d rows but %s image has %d volumes. Regenerate the whole-brain maps and metadata together.', ...
        prefix, height(metadata_table), object_name, n_images);
end
end

% =========================================================================
% Utilities
% =========================================================================
function c = local_to_cell(x)
if isempty(x)
    c = {};
elseif isa(x, 'image_vector')
    c = {x};
elseif ischar(x) || isstring(x)
    c = cellstr(string(x));
else
    c = x;
end
end

function x = local_to_numeric(x)
if isnumeric(x)
    x = double(x);
else
    x = str2double(string(x));
end
end

function context = local_warning_context(user_context, model_name, object_name, prefix)
parts = {};
user_context = char(user_context);
if ~isempty(user_context)
    parts{end + 1} = user_context;
end
if ~isempty(model_name)
    parts{end + 1} = sprintf('model=%s', model_name);
end
parts{end + 1} = sprintf('object=%s', char(object_name));
parts{end + 1} = sprintf('prefix=%s', char(prefix));
context = strjoin(parts, '; ');
end
