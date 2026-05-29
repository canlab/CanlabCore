function avg = hrf_make_average_montage_animations(input_data, varargin)
%HRF_MAKE_AVERAGE_MONTAGE_ANIMATIONS Average whole-brain HRF maps and animate.
%
% avg = hrf_make_average_montage_animations(input_data, ...)
%
% input_data can be:
%   - a table from hrf_collect_wholebrain_outputs
%   - a CSV filename for that table
%   - a study struct from hrf_input_table_to_study with loaded wholebrain maps
%
% The function averages maps across runs within each subject for each
% condition, then averages those subject means across subjects. If OutputDir
% is supplied, subject-level and group-level montage animations are written
% with hrf_make_montage_animation.

p = inputParser;
p.addRequired('input_data', @(x) istable(x) || isstruct(x) || ischar(x) || isstring(x));
p.addParameter('Object', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('Model', '', @(x) ischar(x) || isstring(x));
p.addParameter('Conditions', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('MakeAnimations', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MakeSubjectAnimations', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MakeGroupAnimations', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('ReturnImages', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('GroupStatistic', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('GroupCorrection', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('GroupAlpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('FrameStep', 1, @(x) isscalar(x) && x >= 1);
p.addParameter('FPS', 8, @(x) isscalar(x) && x > 0);
p.addParameter('Threshold', [], @(x) isempty(x) || isscalar(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('SubjectThreshold', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('SubjectCorrection', '', @(x) ischar(x) || isstring(x));
p.addParameter('SubjectAlpha', NaN, @(x) isscalar(x));
p.addParameter('ColorLimits', 'auto', @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && numel(x) == 2));
p.addParameter('ColorMap', 'rdylbu', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('GrayBuffer', false, @(x) islogical(x) || isnumeric(x));
p.parse(input_data, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);

% Resolve subject-level thresholding defaults from the group settings so
% they track unless the user overrides explicitly.
subject_correction = char(opts.SubjectCorrection);
if isempty(subject_correction)
    subject_correction = char(opts.GroupCorrection);
end
subject_alpha = opts.SubjectAlpha;
if isnan(subject_alpha)
    subject_alpha = opts.GroupAlpha;
end
opts.SubjectCorrectionResolved = subject_correction;
opts.SubjectAlphaResolved = subject_alpha;

records = local_records(input_data);
records = local_filter_records(records, char(opts.Model));
if isempty(records)
    error('No whole-brain records were available after filtering.');
end

% Per-prefix cache for the loaded wholebrain struct (and fallback single-file
% loads). containers.Map is a handle class -- mutating it inside the local
% helpers mutates the same instance in the main loop, so each unique source
% NIfTI is read from disk at most once across all conditions, subjects, and
% runs. On UNC-mounted study directories this is the dominant speedup.
record_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

condition_patterns = local_condition_patterns(records, opts.Conditions, record_cache);
output_dir = char(opts.OutputDir);
write_movies = logical(opts.MakeAnimations) && ~isempty(output_dir);
if write_movies && exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end

subject_ids = unique({records.subject}, 'stable');
subject_rows = {};
group_rows = {};

if verbose
    fprintf('\nhrf_make_average_montage_animations\n');
    fprintf('  %d records  |  %d conditions  |  %d subjects  |  model=%s, object=%s\n', ...
        numel(records), numel(condition_patterns), numel(subject_ids), ...
        char(opts.Model), char(opts.Object));
    if write_movies
        fprintf('  output_dir: %s\n', output_dir);
    end
    fprintf('\n');
end
start_time = tic;

avg = struct();
avg.object = char(opts.Object);
avg.model = char(opts.Model);
avg.conditions = condition_patterns(:);
avg.output_dir = output_dir;
avg.subject = struct();
avg.group = struct();
avg.group_t = struct();

for c = 1:numel(condition_patterns)
    condition_pattern = condition_patterns{c};
    subject_images = {};
    subject_metas = {};
    subject_image_subjects = {};
    cond_time = tic;

    if verbose
        fprintf('[%d/%d] condition = %s\n', c, numel(condition_patterns), condition_pattern);
    end

    for s = 1:numel(subject_ids)
        sid = subject_ids{s};
        subject_mask = strcmp({records.subject}, sid);
        subject_records = records(subject_mask(:));
        cache_before = record_cache.Count;
        subj_time = tic;
        [subject_img, subject_meta, n_runs, skipped] = local_subject_average( ...
            subject_records, condition_pattern, opts, record_cache);
        dt = toc(subj_time);
        new_loads = record_cache.Count - cache_before;
        if verbose
            tag = local_load_tag(new_loads, numel(subject_records));
            fprintf('  [%2d/%2d] %-12s %d runs  %s  %5.1fs\n', ...
                s, numel(subject_ids), sid, n_runs, tag, dt);
        end

        if local_missing_image(subject_img)
            local_handle_missing(opts.MissingPolicy, ...
                sprintf('No valid runs for subject %s, condition %s.', sid, condition_pattern));
            continue
        end

        condition_label = local_condition_label(subject_meta, condition_pattern);
        movie_file = '';
        if write_movies && logical(opts.MakeSubjectAnimations)
            movie_file = fullfile(output_dir, sprintf('%s_%s_%s_subject_average.mp4', ...
                local_safe_label(sid), local_safe_label(condition_label), local_object_label(opts.Object)));
            hrf_make_montage_animation(subject_img, movie_file, ...
                'FrameStep', opts.FrameStep, ...
                'FPS', opts.FPS, ...
                'Threshold', opts.Threshold, ...
                'ColorLimits', opts.ColorLimits, ...
                'ColorMap', opts.ColorMap, ...
                'GrayBuffer', opts.GrayBuffer, ...
                'UseSigMask', true, ...
                'TitlePrefix', sprintf('%s %s', sid, condition_label));
        end

        if logical(opts.ReturnImages)
            avg.subject.(matlab.lang.makeValidName(sid)).(matlab.lang.makeValidName(condition_label)) = subject_img;
        end
        subject_rows(end + 1, :) = {sid, condition_label, condition_pattern, char(opts.Object), n_runs, movie_file, numel(skipped)}; %#ok<AGROW>
        subject_images{end + 1, 1} = subject_img; %#ok<AGROW>
        subject_metas{end + 1, 1} = subject_meta; %#ok<AGROW>
        subject_image_subjects{end + 1, 1} = sid; %#ok<AGROW>
    end

    if verbose
        fprintf('  group: averaging across %d subjects ...\n', numel(subject_images));
    end
    group_time = tic;
    [group_img, group_t_img, group_meta, n_subjects] = local_group_average(subject_images, subject_metas, condition_pattern, opts);
    if local_missing_image(group_img)
        local_handle_missing(opts.MissingPolicy, ...
            sprintf('No valid subjects for group condition %s.', condition_pattern));
        continue
    end
    if verbose
        [n_sig_v_preview, ~] = local_sig_counts(group_t_img);
        fprintf('  group: %d subjects, %d significant voxel-frames  %5.1fs\n', ...
            n_subjects, n_sig_v_preview, toc(group_time));
    end

    group_condition_label = local_condition_label(group_meta, condition_pattern);
    group_movie_img = local_group_movie_image(group_img, group_t_img, opts);
    group_movie_file = '';
    if write_movies && logical(opts.MakeGroupAnimations)
        group_movie_file = fullfile(output_dir, sprintf('group_%s_%s_average.mp4', ...
            local_safe_label(group_condition_label), local_group_stat_label(opts)));
        hrf_make_montage_animation(group_movie_img, group_movie_file, ...
            'FrameStep', opts.FrameStep, ...
            'FPS', opts.FPS, ...
            'Threshold', opts.Threshold, ...
            'ColorLimits', opts.ColorLimits, ...
            'ColorMap', opts.ColorMap, ...
            'GrayBuffer', opts.GrayBuffer, ...
            'UseSigMask', true, ...
            'TitlePrefix', sprintf('Group %s', group_condition_label));
    end

    if logical(opts.ReturnImages)
        avg.group.(matlab.lang.makeValidName(group_condition_label)) = group_img;
        avg.group_t.(matlab.lang.makeValidName(group_condition_label)) = group_t_img;
    end
    [n_sig_voxels, n_sig_timepoints] = local_sig_counts(group_t_img);
    group_rows(end + 1, :) = {group_condition_label, condition_pattern, char(opts.Object), ...
        local_group_stat_label(opts), lower(char(opts.GroupCorrection)), opts.GroupAlpha, ...
        n_subjects, group_movie_file, strjoin(subject_image_subjects, ','), ...
        n_sig_voxels, n_sig_timepoints}; %#ok<AGROW>

    if verbose
        fprintf('  condition done in %5.1fs  (cache: %d entries)\n\n', ...
            toc(cond_time), record_cache.Count);
    end
end

if verbose
    fprintf('hrf_make_average_montage_animations complete in %5.1fs  (cache: %d entries)\n\n', ...
        toc(start_time), record_cache.Count);
end

avg.subject_table = local_subject_table(subject_rows);
avg.group_table = local_group_table(group_rows);
end


function tag = local_load_tag(new_loads, n_records)
% Compact indicator: '[cache hit]', '[load 3/3]', or '[load 1/3, hit 2/3]'.
if new_loads == 0
    tag = '[cache hit]';
elseif new_loads == n_records
    tag = sprintf('[load %d/%d]', new_loads, n_records);
else
    tag = sprintf('[load %d/%d, hit %d/%d]', new_loads, n_records, n_records - new_loads, n_records);
end
end

function tf = local_missing_image(img)
if isobject(img) && isprop(img, 'dat')
    tf = isempty(img.dat);
else
    tf = isempty(img);
end
end

function records = local_records(input_data)
if istable(input_data) || ischar(input_data) || isstring(input_data)
    records = local_table_records(local_read_table(input_data));
    return
end

records = local_study_records(input_data);
if isempty(records) && isfield(input_data, 'second_level_inputs')
    records = local_table_records(input_data.second_level_inputs);
end
end

function T = local_read_table(input_data)
if istable(input_data)
    T = input_data;
else
    T = readtable(char(input_data), 'TextType', 'string');
end
end

function records = local_table_records(T)
records = repmat(local_empty_record(), 0, 1);
for i = 1:height(T)
    rec = local_empty_record();
    rec.subject = local_table_value(T, i, 'subject');
    rec.run_label = local_table_value(T, i, 'run_label');
    rec.model = lower(local_table_value(T, i, 'model'));
    rec.prefix = local_table_value(T, i, 'prefix');
    rec.beta_file = local_table_value(T, i, 'beta_file');
    rec.t_file = local_table_value(T, i, 't_file');
    rec.metadata_file = local_table_value(T, i, 'metadata_file');
    if isempty(rec.subject)
        rec.subject = sprintf('sub-%03d', i);
    end
    if isempty(rec.run_label)
        rec.run_label = sprintf('run-%03d', i);
    end
    records(end + 1, 1) = rec; %#ok<AGROW>
end
end

function records = local_study_records(study)
records = repmat(local_empty_record(), 0, 1);
if ~isfield(study, 'results')
    return
end

for i = 1:numel(study.results)
    r = study.results{i};
    if isempty(r)
        continue
    end
    rec = local_empty_record();
    rec.subject = local_subject_id(study, r, i);
    rec.run_label = local_run_label(study, r, i, rec.subject);
    if isfield(r, 'wholebrain') && ~isempty(r.wholebrain)
        rec.wholebrain = r.wholebrain;
        rec.model = local_model_from_metadata(r.wholebrain);
    end
    if isfield(r, 'input_prefix')
        rec.prefix = char(r.input_prefix);
    end
    if isempty(rec.model) && isfield(study, 'second_level_inputs') && height(study.second_level_inputs) >= i
        rec.model = lower(local_table_value(study.second_level_inputs, i, 'model'));
    end
    records(end + 1, 1) = rec; %#ok<AGROW>
end
end

function rec = local_empty_record()
rec = struct('subject', '', 'run_label', '', 'model', '', 'prefix', '', ...
    'beta_file', '', 't_file', '', 'metadata_file', '', 'wholebrain', []);
end

function records = local_filter_records(records, model_name)
model_name = lower(strtrim(model_name));
if isempty(model_name)
    return
end
keep = false(size(records));
for i = 1:numel(records)
    keep(i) = isempty(records(i).model) || strcmpi(records(i).model, model_name);
end
records = records(keep);
end

function patterns = local_condition_patterns(records, requested, cache)
if ~isempty(requested)
    patterns = cellstr(string(requested));
    patterns = patterns(:)';
    return
end

patterns = {};
for i = 1:numel(records)
    try
        [~, metadata_table] = local_record_object(records(i), 'beta', cache);
    catch
        continue
    end
    if ~isempty(metadata_table) && any(strcmp('condition', metadata_table.Properties.VariableNames))
        patterns = [patterns; cellstr(string(metadata_table.condition))]; %#ok<AGROW>
    end
end
patterns = unique(patterns, 'stable');
if isempty(patterns)
    error('Could not infer conditions. Supply the Conditions option.');
end
patterns = patterns(:)';
end

function [subject_img, subject_meta, n_runs, skipped] = local_subject_average(records, condition_pattern, opts, cache)
run_dat = {};
template_img = [];
subject_meta = table();
skipped = {};

for i = 1:numel(records)
    try
        [obj, metadata_table] = local_record_object(records(i), char(opts.Object), cache);
        [Y, meta] = local_condition_matrix(obj, metadata_table, condition_pattern);
        if isempty(template_img)
            template_img = obj;
            subject_meta = meta;
        else
            local_assert_same_time(subject_meta, meta);
        end
        run_dat{end + 1, 1} = Y; %#ok<AGROW>
    catch err
        skipped{end + 1, 1} = sprintf('%s: %s', records(i).run_label, err.message); %#ok<AGROW>
        if strcmpi(char(opts.MissingPolicy), 'error')
            rethrow(err);
        elseif strcmpi(char(opts.MissingPolicy), 'warn')
            warning('hrf_make_average_montage_animations:SkippingRun', ...
                'Skipping %s %s: %s', records(i).subject, records(i).run_label, err.message);
        end
    end
end

n_runs = numel(run_dat);
if n_runs == 0
    subject_img = [];
    return
end

stack = local_stack(run_dat);
[mean_dat, sem_dat, n_dat] = local_mean_sem_stack(stack);
subject_img = local_make_average_image(template_img, mean_dat, sem_dat, subject_meta, ...
    sprintf('Subject average, n runs = %d', n_runs), n_runs);

% Optionally populate subject_img.sig from a per-subject t-stat (mean/sem)
% so the animation's UseSigMask actually masks something. n_runs < 2
% means no SE is available -- skip.
if logical(opts.SubjectThreshold) && n_runs >= 2 && isa(subject_img, 'statistic_image')
    correction = opts.SubjectCorrectionResolved;
    alpha = opts.SubjectAlphaResolved;
    subject_t_img = local_make_t_image(template_img, mean_dat, sem_dat, n_dat, ...
        subject_meta, correction, alpha, ...
        sprintf('Subject one-sample t maps across runs, n runs = %d', n_runs));
    if isa(subject_t_img, 'statistic_image') && ~isempty(subject_t_img.sig)
        subject_img.sig = subject_t_img.sig;
        subject_img.p = subject_t_img.p;
        subject_img.p_type = subject_t_img.p_type;
        subject_img.ste = sem_dat;
    end
end
end

function [group_img, group_t_img, group_meta, n_subjects] = local_group_average(subject_images, subject_metas, condition_pattern, opts)
group_img = [];
group_t_img = [];
group_meta = table();
n_subjects = numel(subject_images);
if n_subjects == 0
    return
end

subject_dat = cell(n_subjects, 1);
template = subject_images{1};
for i = 1:n_subjects
    subject_dat{i} = subject_images{i}.dat;
end
stack = local_stack(subject_dat);
[mean_dat, sem_dat, n_dat] = local_mean_sem_stack(stack);
if ~isempty(subject_metas) && ~isempty(subject_metas{1})
    group_meta = subject_metas{1};
else
    group_meta = local_metadata_from_image(template, condition_pattern);
end
group_img = local_make_average_image(template, mean_dat, sem_dat, group_meta, ...
    sprintf('Group average, n subjects = %d', n_subjects), n_subjects);
group_t_img = local_make_t_image(template, mean_dat, sem_dat, n_dat, group_meta, ...
    opts.GroupCorrection, opts.GroupAlpha, ...
    sprintf('Group one-sample t maps across subject averages, n subjects = %d', max(n_dat(:))));
end

function [obj, metadata_table] = local_record_object(record, object_name, cache)
% Cache strategy:
%   - record.wholebrain in memory      -> no disk read, no cache needed
%   - record.prefix given              -> cache the loaded wholebrain struct
%                                         under 'wb:<prefix>'; serves all
%                                         object_name requests for free
%   - direct single-file fallback      -> cache the (obj, metadata) tuple
%                                         under 'f:<file>'
% containers.Map is a handle class, so mutations propagate to all callers.

if ~isempty(record.wholebrain)
    [obj, metadata_table] = local_object_from_wholebrain(record.wholebrain, object_name);
    return
end

if ~isempty(record.prefix)
    wb_key = ['wb:' record.prefix];
    if isKey(cache, wb_key)
        wb = cache(wb_key);
        [obj, metadata_table] = local_object_from_wholebrain(wb, object_name);
        return
    end
    try
        wb = hrf_load_wholebrain_stats(record.prefix, 'NoVerbose', true);
        cache(wb_key) = wb;
        [obj, metadata_table] = local_object_from_wholebrain(wb, object_name);
        return
    catch
    end
end

[file_name, type_label] = local_object_file(record, object_name);
if isempty(file_name) || exist(file_name, 'file') ~= 2
    error('Missing %s image file.', type_label);
end
f_key = ['f:' file_name];
if isKey(cache, f_key)
    cached = cache(f_key);
    obj = cached.obj;
    metadata_table = cached.metadata;
    return
end
obj = statistic_image(fmri_data(file_name, 'noverbose'));
obj.type = type_label;
metadata_table = local_read_metadata(record.metadata_file);
obj = local_attach_metadata(obj, metadata_table);
cache(f_key) = struct('obj', obj, 'metadata', metadata_table);
end

function [obj, metadata_table] = local_object_from_wholebrain(wholebrain, object_name)
switch lower(char(object_name))
    case {'beta', 'b'}
        obj = wholebrain.b;
    case {'t', 'tmap', 'tmaps'}
        obj = wholebrain.t;
    otherwise
        error('Unknown Object: %s. Use ''beta'' or ''t''.', char(object_name));
end
if isfield(wholebrain, 'metadata_table')
    metadata_table = wholebrain.metadata_table;
else
    metadata_table = table();
end
end

function [file_name, type_label] = local_object_file(record, object_name)
switch lower(char(object_name))
    case {'beta', 'b'}
        file_name = record.beta_file;
        type_label = 'beta';
    case {'t', 'tmap', 'tmaps'}
        file_name = record.t_file;
        type_label = 't';
    otherwise
        error('Unknown Object: %s. Use ''beta'' or ''t''.', char(object_name));
end
end

function metadata_table = local_read_metadata(metadata_file)
if ~isempty(metadata_file) && exist(metadata_file, 'file') == 2
    metadata_table = readtable(metadata_file, 'TextType', 'string');
else
    metadata_table = table();
end
end

function obj = local_attach_metadata(obj, metadata_table)
if isempty(metadata_table) || height(metadata_table) ~= size(obj.dat, 2)
    return
end
if any(strcmp('image_label', metadata_table.Properties.VariableNames))
    obj.image_labels = cellstr(string(metadata_table.image_label));
end
if any(strcmp('N', metadata_table.Properties.VariableNames))
    obj.N = metadata_table.N(1);
end
if any(strcmp('dfe', metadata_table.Properties.VariableNames))
    obj.dfe = metadata_table.dfe(1);
end
if isempty(obj.sig)
    obj.sig = true(size(obj.dat));
end
end

function [Y, meta_out] = local_condition_matrix(obj, metadata_table, condition_pattern)
if isempty(metadata_table) || ~any(strcmp('condition', metadata_table.Properties.VariableNames))
    error('Condition selection needs metadata_table.condition.');
end

available = unique(cellstr(string(metadata_table.condition)), 'stable');
specs = hrf_resolve_condition_patterns(available, condition_pattern, 'DefaultMode', 'first');
spec = specs(1);
cond = cellstr(string(metadata_table.condition));
keep = ismember(cond, spec.matched_conditions);
if ~any(keep)
    error('Condition pattern "%s" matched no image volumes.', condition_pattern);
end

if any(strcmp('lag_index', metadata_table.Properties.VariableNames))
    [Y, meta_out] = local_condition_lag_average(obj, metadata_table, keep, spec);
else
    wh = find(keep);
    Y = obj.dat(:, wh);
    meta_out = metadata_table(wh, :);
end
end

function [Y, meta_out] = local_condition_lag_average(obj, metadata_table, keep, spec)
lags = unique(metadata_table.lag_index(keep), 'stable');
Y = nan(size(obj.dat, 1), numel(lags));
rows = cell(numel(lags), width(metadata_table));
var_names = metadata_table.Properties.VariableNames;

for lag_i = 1:numel(lags)
    wh = keep & metadata_table.lag_index == lags(lag_i);
    Y(:, lag_i) = local_mean_omitnan(obj.dat(:, wh), 2);
    row = metadata_table(find(wh, 1), :);
    if any(strcmp('condition', var_names))
        row = local_set_row_value(row, 'condition', spec.display_label);
    end
    if any(strcmp('image_label', var_names))
        row = local_set_row_value(row, 'image_label', local_lag_label(spec.label, row, lags(lag_i)));
    end
    rows(lag_i, :) = table2cell(row);
end

function row = local_set_row_value(row, varname, value)
current = row.(varname);
if iscell(current)
    row.(varname) = {char(value)};
elseif isstring(current)
    row.(varname) = string(value);
elseif iscategorical(current)
    row.(varname) = categorical({char(value)});
else
    row.(varname) = string(value);
end
end
meta_out = cell2table(rows, 'VariableNames', var_names);
if any(strcmp('volume_index', var_names))
    meta_out.volume_index = (1:height(meta_out))';
end
if any(strcmp('condition_index', var_names))
    meta_out.condition_index = ones(height(meta_out), 1);
end
end

function label = local_lag_label(condition_label, row, lag_index)
safe_condition = local_safe_label(condition_label);
if any(strcmp('lag_seconds', row.Properties.VariableNames))
    label = sprintf('%s_lag%03d_%0.3gs', safe_condition, lag_index, row.lag_seconds(1));
else
    label = sprintf('%s_lag%03d', safe_condition, lag_index);
end
end

function local_assert_same_time(a, b)
if height(a) ~= height(b)
    error('Condition time-bin count mismatch.');
end
if any(strcmp('lag_index', a.Properties.VariableNames)) && any(strcmp('lag_index', b.Properties.VariableNames)) && ...
        ~isequal(a.lag_index, b.lag_index)
    error('Condition lag_index mismatch.');
end
if any(strcmp('lag_seconds', a.Properties.VariableNames)) && any(strcmp('lag_seconds', b.Properties.VariableNames)) && ...
        any(abs(a.lag_seconds - b.lag_seconds) > max(eps(a.lag_seconds)))
    error('Condition lag_seconds mismatch.');
end
end

function stack = local_stack(dat_cells)
n = numel(dat_cells);
sz = size(dat_cells{1});
stack = nan([sz n]);
for i = 1:n
    if ~isequal(size(dat_cells{i}), sz)
        error('Image data size mismatch while averaging.');
    end
    stack(:, :, i) = dat_cells{i};
end
end

function [m, se, n] = local_mean_sem_stack(X)
valid = ~isnan(X);
n = sum(valid, 3);
X0 = X;
X0(~valid) = 0;
m = sum(X0, 3) ./ n;
m(n == 0) = NaN;

centered = X - repmat(m, 1, 1, size(X, 3));
centered(~valid) = 0;
sd = sqrt(sum(centered .^ 2, 3) ./ max(n - 1, 1));
se = sd ./ sqrt(n);
se(n == 0) = NaN;
end

function t_img = local_make_t_image(template_img, mean_dat, sem_dat, n_dat, metadata_table, correction, alpha, description)
% One-sample t map with corrected sig mask. Shared between subject- and
% group-level paths (n_dat is "samples per voxel-volume" -- subjects at
% group level, runs at subject level).
t_dat = mean_dat ./ sem_dat;
t_dat(~isfinite(t_dat) | n_dat < 2) = NaN;
dfe = max(n_dat - 1, 1);
p_dat = 2 * (1 - tcdf(abs(t_dat), dfe));
p_dat(n_dat < 2 | ~isfinite(t_dat)) = NaN;

t_img = template_img;
t_img.dat = t_dat;
if isa(t_img, 'statistic_image')
    t_img.type = 'T';
    t_img.p = p_dat;
    t_img.p_type = sprintf('Two-tailed one-sample t-test, %s q < %.3g', ...
        upper(char(correction)), alpha);
    t_img.ste = sem_dat;
    t_img.sig = local_corrected_sig(p_dat, alpha, correction);
    t_img.dat_descrip = char(description);
end
if isprop(t_img, 'image_labels') && ~isempty(metadata_table) && any(strcmp('image_label', metadata_table.Properties.VariableNames))
    t_img.image_labels = cellstr(string(metadata_table.image_label));
end
if isprop(t_img, 'removed_images')
    t_img.removed_images = false(1, size(t_dat, 2));
end
if isprop(t_img, 'N')
    t_img.N = max(n_dat(:));
end
if isprop(t_img, 'dfe')
    t_img.dfe = max(max(n_dat(:)) - 1, 1);
end
end

function sig = local_corrected_sig(p_dat, alpha, correction)
sig = false(size(p_dat));
correction = lower(strtrim(char(correction)));
switch correction
    case {'none', 'unc', 'uncorrected'}
        sig = p_dat < alpha;
    case 'fdr'
        for i = 1:size(p_dat, 2)
            p = p_dat(:, i);
            valid = isfinite(p);
            if ~any(valid)
                continue
            end
            thr = FDR(p(valid), alpha);
            if isempty(thr) || ~isfinite(thr)
                continue
            end
            sig(valid, i) = p(valid) <= thr;
        end
    otherwise
        error('Unknown GroupCorrection: %s. Use ''none'' or ''fdr''.', char(correction));
end
sig(isnan(p_dat)) = false;
end

function img = local_group_movie_image(group_img, group_t_img, opts)
switch lower(strtrim(char(opts.GroupStatistic)))
    case {'mean', 'average', 'beta'}
        img = group_img;
    case {'t', 'tmap', 'tstat'}
        img = group_t_img;
    otherwise
        error('Unknown GroupStatistic: %s. Use ''mean'' or ''t''.', char(opts.GroupStatistic));
end
end

function label = local_group_stat_label(opts)
switch lower(strtrim(char(opts.GroupStatistic)))
    case {'mean', 'average', 'beta'}
        label = local_object_label(opts.Object);
    case {'t', 'tmap', 'tstat'}
        correction = strtrim(char(opts.GroupCorrection));
        if strcmpi(correction, 'fdr')
            label = sprintf('t_fdr_q%0.3g', opts.GroupAlpha);
        else
            label = 't';
        end
        label = local_safe_label(label);
    otherwise
        label = local_safe_label(opts.GroupStatistic);
end
end

function [n_sig_voxels, n_sig_timepoints] = local_sig_counts(group_t_img)
n_sig_voxels = 0;
n_sig_timepoints = 0;
if isa(group_t_img, 'statistic_image') && ~isempty(group_t_img.sig)
    n_sig_voxels = sum(group_t_img.sig(:));
    n_sig_timepoints = sum(any(group_t_img.sig, 1));
end
end

function img = local_make_average_image(template_img, dat, ste, metadata_table, description, n_obs)
img = template_img;
img.dat = dat;
if isa(img, 'statistic_image')
    img.p = [];
    img.sig = true(size(dat));
    img.ste = ste;
    img.dat_descrip = description;
end
if isprop(img, 'image_labels') && ~isempty(metadata_table) && any(strcmp('image_label', metadata_table.Properties.VariableNames))
    img.image_labels = cellstr(string(metadata_table.image_label));
end
if isprop(img, 'removed_images')
    img.removed_images = false(1, size(dat, 2));
end
if isprop(img, 'N')
    img.N = n_obs;
end
if isprop(img, 'dfe')
    img.dfe = [];
end
end

function metadata_table = local_metadata_from_image(img, condition_pattern)
metadata_table = table();
if isprop(img, 'image_labels') && ~isempty(img.image_labels)
    labels = cellstr(string(img.image_labels(:)));
    n = numel(labels);
    metadata_table = table((1:n)', repmat(string(condition_pattern), n, 1), (1:n)', labels(:), ...
        'VariableNames', {'volume_index', 'condition', 'lag_index', 'image_label'});
end
end

function label = local_condition_label(metadata_table, fallback)
label = char(fallback);
if ~isempty(metadata_table) && any(strcmp('condition', metadata_table.Properties.VariableNames)) && height(metadata_table) > 0
    label = char(string(metadata_table.condition(1)));
end
end

function T = local_subject_table(rows)
vars = {'subject', 'condition', 'condition_pattern', 'object', 'n_runs', 'movie_file', 'n_skipped_runs'};
if isempty(rows)
    T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
else
    T = cell2table(rows, 'VariableNames', vars);
end
end

function T = local_group_table(rows)
vars = {'condition', 'condition_pattern', 'object', 'group_statistic', 'group_correction', ...
    'group_alpha', 'n_subjects', 'movie_file', 'subjects', 'n_sig_voxels', 'n_sig_timepoints'};
if isempty(rows)
    T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
else
    T = cell2table(rows, 'VariableNames', vars);
end
end

function local_handle_missing(missing_policy, message_text)
switch lower(char(missing_policy))
    case 'error'
        error(message_text);
    case 'warn'
        warning('hrf_make_average_montage_animations:MissingAverage', '%s', message_text);
    case 'silent'
    otherwise
        error('Unknown MissingPolicy: %s. Use ''warn'', ''silent'', or ''error''.', char(missing_policy));
end
end

function subject_id = local_subject_id(study, r, idx)
if isfield(r, 'subject_id') && ~isempty(r.subject_id)
    subject_id = char(r.subject_id);
elseif isfield(study, 'subject_ids') && numel(study.subject_ids) >= idx
    subject_id = char(study.subject_ids{idx});
else
    subject_id = sprintf('sub-%03d', idx);
end
end

function run_label = local_run_label(study, r, idx, subject_id)
if isfield(r, 'run_label') && ~isempty(r.run_label)
    run_label = char(r.run_label);
elseif isfield(study, 'run_labels') && numel(study.run_labels) >= idx && ~isempty(study.run_labels{idx})
    run_label = char(study.run_labels{idx});
else
    run_label = sprintf('%s_run%03d', subject_id, idx);
end
end

function model_name = local_model_from_metadata(wholebrain)
model_name = '';
if isfield(wholebrain, 'metadata_table') && ~isempty(wholebrain.metadata_table) && ...
        any(strcmp('mode', wholebrain.metadata_table.Properties.VariableNames)) && height(wholebrain.metadata_table) > 0
    model_name = lower(char(string(wholebrain.metadata_table.mode(1))));
end
end

function value = local_table_value(T, row, varname)
if ~any(strcmp(varname, T.Properties.VariableNames))
    value = '';
    return
end
value = T.(varname)(row);
if iscell(value)
    value = value{1};
end
value = string(value);
if ismissing(value) || strlength(value) == 0
    value = '';
else
    value = char(value);
end
end

function s = local_safe_label(s)
s = regexprep(char(s), '[*?]', 'all');
s = matlab.lang.makeValidName(s);
end

function s = local_object_label(object_name)
switch lower(char(object_name))
    case {'beta', 'b'}
        s = 'beta';
    case {'t', 'tmap', 'tmaps'}
        s = 't';
    otherwise
        s = local_safe_label(object_name);
end
end

function m = local_mean_omitnan(X, dim)
if nargin < 2
    dim = 1;
end
valid = ~isnan(X);
den = sum(valid, dim);
X(~valid) = 0;
m = sum(X, dim) ./ den;
m(den == 0) = NaN;
end
