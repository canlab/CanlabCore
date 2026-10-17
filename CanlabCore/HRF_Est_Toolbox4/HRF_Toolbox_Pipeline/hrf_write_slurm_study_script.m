function paths = hrf_write_slurm_study_script(fmri_files, events_files, subject_ids, varargin)
%HRF_WRITE_SLURM_STUDY_SCRIPT Write a SLURM array job for HRF study fitting.
%
% paths = hrf_write_slurm_study_script(fmri_files, events_files, subject_ids, ...)
%
% The generated job runs one fMRI/events pair per SLURM array task. Each
% task calls run_hrf_pipeline, writes 4D whole-brain beta/T statistic_image
% files, and optionally applies all requested signatures/imagesets to those
% 4D maps.

if ischar(fmri_files) || isstring(fmri_files), fmri_files = cellstr(fmri_files); end
if ischar(events_files) || isstring(events_files), events_files = cellstr(events_files); end
if nargin < 3 || isempty(subject_ids)
    subject_ids = arrayfun(@(i) sprintf('sub-%03d', i), 1:numel(fmri_files), 'UniformOutput', false);
end
if isstring(subject_ids), subject_ids = cellstr(subject_ids); end

n = numel(fmri_files);
if numel(events_files) ~= n
    error('fmri_files and events_files must contain the same number of files.');
end
if numel(subject_ids) ~= n
    error('subject_ids must be empty or contain one id per fMRI file.');
end

p = inputParser;
p.addParameter('OutputDir', fullfile(pwd, 'hrf_outputs'), @(x) ischar(x) || isstring(x));
p.addParameter('ScriptFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('ManifestFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('WorkerFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConfigMat', '', @(x) ischar(x) || isstring(x));
p.addParameter('RunLabels', {}, @(x) isempty(x) || ischar(x) || iscell(x) || isstring(x));
p.addParameter('CanlabRoot', local_default_canlab_root(), @(x) ischar(x) || isstring(x));
p.addParameter('ExtraMatlabPaths', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('PipelineArgs', {}, @(x) iscell(x));
% Per-subject estimated SPM.mat paths for exact Tier B g*K*W matching (one
% per fMRI file; '' = none for that row). SPMRuns picks the run within a
% multi-run SPM (default 1 each). See hrf_fit_wholebrain_stats.
p.addParameter('SPMFiles', {}, @(x) isempty(x) || ischar(x) || iscell(x) || isstring(x));
p.addParameter('SPMRuns', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('SignatureSets', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ImageSets', {}, @(x) ischar(x) || iscell(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('AtlasObj', [], @(x) isempty(x) || isa(x, 'atlas') || isa(x, 'image_vector'));
p.addParameter('AtlasName', '', @(x) ischar(x) || isstring(x));
p.addParameter('Regions', {}, @(x) iscell(x) || isstring(x) || ischar(x));
% Multi-atlas: score several atlases in one job. Give each a DISTINCT name in
% AtlasNames so their region columns do not collide. AtlasRegions is a cell of
% per-atlas region lists ({} = all regions for that atlas).
p.addParameter('AtlasObjs', {}, @(x) isempty(x) || iscell(x) || isa(x, 'atlas') || isa(x, 'image_vector'));
p.addParameter('AtlasNames', {}, @(x) isempty(x) || iscell(x) || ischar(x) || isstring(x));
p.addParameter('AtlasRegions', {}, @(x) isempty(x) || iscell(x));
p.addParameter('Normalize', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('ScoreObjects', {'beta'}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.addParameter('MakeAnimations', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('AnimationObject', 't', @(x) ischar(x) || isstring(x));
p.addParameter('AnimationCondition', '', @(x) ischar(x) || isstring(x));
p.addParameter('AnimationFrameRate', 4, @(x) isscalar(x) && x > 0);
p.addParameter('JobName', 'hrf_study', @(x) ischar(x) || isstring(x));
p.addParameter('Time', '24:00:00', @(x) ischar(x) || isstring(x));
p.addParameter('Mem', '32G', @(x) ischar(x) || isstring(x));
p.addParameter('CpusPerTask', 1, @(x) isscalar(x) && x >= 1);
p.addParameter('Partition', '', @(x) ischar(x) || isstring(x));
p.addParameter('Account', '', @(x) ischar(x) || isstring(x));
p.addParameter('LogDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('ModuleLoad', '', @(x) ischar(x) || isstring(x));
p.addParameter('MatlabCommand', 'matlab', @(x) ischar(x) || isstring(x));
p.addParameter('ExtraSlurmDirectives', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.parse(varargin{:});
opts = p.Results;

output_dir = char(opts.OutputDir);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

log_dir = char(opts.LogDir);
if isempty(log_dir)
    log_dir = fullfile(output_dir, 'logs');
end
if ~exist(log_dir, 'dir'), mkdir(log_dir); end

paths = local_default_paths(output_dir, opts);
local_ensure_parent(paths.script_file);
local_ensure_parent(paths.manifest_file);
local_ensure_parent(paths.worker_file);
local_ensure_parent(paths.config_mat);

spm_files = local_normalize_spm_files(opts.SPMFiles, n);
spm_runs = local_normalize_spm_runs(opts.SPMRuns, n);
manifest = local_manifest_table(fmri_files, events_files, subject_ids, output_dir, opts.RunLabels, spm_files, spm_runs);
writetable(manifest, paths.manifest_file);

config = struct();
config.canlab_root = char(opts.CanlabRoot);
config.extra_matlab_paths = local_to_cell(opts.ExtraMatlabPaths);
config.subject_ids = subject_ids(:);
config.run_labels = manifest.run_label(:);
config.fmri_files = fmri_files(:);
config.events_files = events_files(:);
config.output_prefixes = manifest.output_prefix(:);
config.output_mats = manifest.output_mat(:);
config.spm_files = spm_files(:);
config.spm_runs = spm_runs(:);
[pipeline_args, pipeline_note] = local_normalize_pipeline_args(opts.PipelineArgs);
config.pipeline_args = pipeline_args;
config.pipeline_note = pipeline_note;
config.signature_sets = local_to_cell(opts.SignatureSets);
config.image_sets = local_to_cell(opts.ImageSets);
config.atlas_obj = opts.AtlasObj;
config.atlas_name = char(opts.AtlasName);
config.regions = local_to_cell(opts.Regions);
config.atlas_objs = local_to_cell(opts.AtlasObjs);
config.atlas_names = local_to_cell(opts.AtlasNames);
config.atlas_regions = opts.AtlasRegions;
config.normalize = char(opts.Normalize);
config.score_objects = local_to_cell(opts.ScoreObjects);
config.similarity_metric = char(opts.SimilarityMetric);
config.make_animations = logical(opts.MakeAnimations);
config.animation_object = char(opts.AnimationObject);
config.animation_condition = char(opts.AnimationCondition);
config.animation_frame_rate = opts.AnimationFrameRate;
save(paths.config_mat, 'config');

local_write_worker(paths.worker_file, paths.manifest_file, paths.config_mat, ...
    char(opts.CanlabRoot), local_to_cell(opts.ExtraMatlabPaths));
local_write_sbatch(paths.script_file, paths.worker_file, n, log_dir, opts);

paths.message = sprintf(['Successfully wrote HRF SLURM study files for %d task(s):\n' ...
    '  sbatch:   %s\n  worker:   %s\n  manifest: %s\n  config:   %s'], ...
    n, paths.script_file, paths.worker_file, paths.manifest_file, paths.config_mat);
if ~isempty(pipeline_note)
    paths.message = sprintf('%s\n%s', paths.message, pipeline_note);
end
fprintf('%s\n', paths.message);
end

function [pipeline_args, note] = local_normalize_pipeline_args(pipeline_args)
note = '';
if isempty(pipeline_args), return; end

models = local_arg_value(pipeline_args, 'Models');
has_wholebrain_mode = local_has_arg(pipeline_args, 'WholeBrainMode');
if ~isempty(models)
    model_names = lower(cellstr(string(models)));
    if ~has_wholebrain_mode
        pipeline_args = [pipeline_args, {'WholeBrainMode', 'auto'}];
        note = sprintf(['Note: PipelineArgs did not include WholeBrainMode, so the SLURM worker will use ' ...
            'WholeBrainMode=''auto'' and write one 4D whole-brain map set for each requested supported linear model.']);
    end

    nonlinear = intersect(model_names, {'logit', 'nlgamma'}, 'stable');
    if ~isempty(nonlinear)
        extra = sprintf(['Note: Models={%s} affects the 1D extracted-signal fits saved in the MAT results, ' ...
            'but nonlinear logit/nlgamma whole-brain maps are skipped by the fast 4D writer.'], ...
            strjoin(nonlinear, ', '));
        if isempty(note), note = extra; else, note = sprintf('%s\n%s', note, extra); end
    end
end
end

function tf = local_has_arg(args, name)
tf = false;
for i = 1:2:numel(args)
    if ischar(args{i}) || isstring(args{i})
        if strcmpi(char(args{i}), name)
            tf = true;
            return
        end
    end
end
end

function value = local_arg_value(args, name)
value = [];
for i = 1:2:numel(args)
    if ischar(args{i}) || isstring(args{i})
        if strcmpi(char(args{i}), name) && i < numel(args)
            value = args{i + 1};
            return
        end
    end
end
end

function paths = local_default_paths(output_dir, opts)
script_file = char(opts.ScriptFile);
if isempty(script_file), script_file = fullfile(output_dir, 'run_hrf_study.sbatch'); end

manifest_file = char(opts.ManifestFile);
if isempty(manifest_file), manifest_file = fullfile(output_dir, 'hrf_study_manifest.csv'); end

worker_file = char(opts.WorkerFile);
if isempty(worker_file), worker_file = fullfile(output_dir, 'run_hrf_study_worker.m'); end

config_mat = char(opts.ConfigMat);
if isempty(config_mat), config_mat = fullfile(output_dir, 'hrf_study_config.mat'); end

paths = struct('script_file', script_file, 'manifest_file', manifest_file, ...
    'worker_file', worker_file, 'config_mat', config_mat);
end

function T = local_manifest_table(fmri_files, events_files, subject_ids, output_dir, run_labels_input, spm_files, spm_runs)
n = numel(fmri_files);
subject = cell(n, 1);
run_label = cell(n, 1);
fmri_file = cell(n, 1);
events_file = cell(n, 1);
output_prefix = cell(n, 1);
output_mat = cell(n, 1);
spm_file = cell(n, 1);
spm_run = zeros(n, 1);
run_labels = local_run_labels(fmri_files, subject_ids, run_labels_input);

for i = 1:n
    subject_label = local_file_label(subject_ids{i});
    run_label{i} = run_labels{i};
    output_label = local_file_label([subject_label '_' run_labels{i}]);
    subject{i} = char(subject_ids{i});
    fmri_file{i} = char(fmri_files{i});
    events_file{i} = char(events_files{i});
    output_prefix{i} = fullfile(output_dir, [output_label '_hrf']);
    output_mat{i} = fullfile(output_dir, [output_label '_hrf_results.mat']);
    spm_file{i} = spm_files{i};
    spm_run(i) = spm_runs(i);
end

T = table((1:n)', subject, run_label, fmri_file, events_file, output_prefix, output_mat, spm_file, spm_run, ...
    'VariableNames', {'index', 'subject', 'run_label', 'fmri_file', 'events_file', 'output_prefix', 'output_mat', 'spm_file', 'spm_run'});
end

function spm_files = local_normalize_spm_files(spm_files, n)
% Per-subject SPM input -> n-element cellstr ('' = none).
if isempty(spm_files)
    spm_files = repmat({''}, 1, n);
    return
end
if ischar(spm_files) || isstring(spm_files)
    spm_files = cellstr(string(spm_files));
end
if isscalar(spm_files)
    spm_files = repmat(spm_files(:)', 1, n);
end
spm_files = cellfun(@(x) char(string(x)), spm_files, 'UniformOutput', false);
if numel(spm_files) ~= n
    error('SPMFiles must be empty, scalar, or contain one entry per fMRI file (got %d, need %d).', numel(spm_files), n);
end
spm_files = reshape(spm_files, 1, n);
end

function spm_runs = local_normalize_spm_runs(spm_runs, n)
if isempty(spm_runs)
    spm_runs = ones(1, n);
elseif isscalar(spm_runs)
    spm_runs = repmat(spm_runs, 1, n);
end
if numel(spm_runs) ~= n
    error('SPMRuns must be empty, scalar, or contain one value per fMRI file.');
end
spm_runs = reshape(spm_runs, 1, n);
end

function run_labels = local_run_labels(fmri_files, subject_ids, run_labels_input)
n = numel(fmri_files);
if ~isempty(run_labels_input)
    run_labels = cellstr(string(run_labels_input));
    if numel(run_labels) ~= n
        error('RunLabels must contain one label per fMRI file.');
    end
    run_labels = cellfun(@local_file_label, run_labels, 'UniformOutput', false);
else
    run_labels = cell(n, 1);
    for i = 1:n
        run_labels{i} = local_run_label_from_filename(fmri_files{i}, i);
    end
end

seen = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:n
    key = sprintf('%s__%s', char(subject_ids{i}), run_labels{i});
    if isKey(seen, key)
        seen(key) = seen(key) + 1;
        run_labels{i} = sprintf('%s_dup%02d', run_labels{i}, seen(key));
    else
        seen(key) = 1;
    end
end
end

function label = local_run_label_from_filename(fmri_file, idx)
[~, name, ext] = fileparts(char(fmri_file));
if strcmpi(ext, '.gz')
    [~, name] = fileparts(name);
end
parts = {};
tokens = {'ses', 'task', 'run', 'acq', 'desc'};
for i = 1:numel(tokens)
    tok = regexp(name, [tokens{i} '-[^_]+'], 'match', 'once');
    if ~isempty(tok)
        parts{end + 1} = tok; %#ok<AGROW>
    end
end
if isempty(parts)
    label = sprintf('run-%03d', idx);
else
    label = strjoin(parts, '_');
end
label = local_file_label(label);
end

function local_write_worker(worker_file, manifest_file, config_mat, canlab_root, extra_paths)
fid = fopen(worker_file, 'w');
if fid == -1, error('Could not write worker file: %s', worker_file); end
c = onCleanup(@() fclose(fid));

fprintf(fid, '%% Auto-generated by hrf_write_slurm_study_script.m\n');
% addpath MUST come before load(config_mat). load() reconstructs classdef
% objects in the saved struct (e.g., the atlas object passed via AtlasObj);
% if the @atlas class folder isn't on path at load-time, the object either
% silently strips to a plain struct or never makes it onto config.atlas_obj
% at all. canlab_root and any ExtraMatlabPaths are known at generation
% time, so they're written as literals.
fprintf(fid, 'addpath(genpath(strrep(''%s'', ''\\'', filesep)));\n', local_matlab_string(canlab_root));
for i = 1:numel(extra_paths)
    p_str = char(string(extra_paths{i}));
    if ~isempty(p_str)
        fprintf(fid, 'addpath(genpath(strrep(''%s'', ''\\'', filesep)));\n', local_matlab_string(p_str));
    end
end
fprintf(fid, 'manifest_file = strrep(''%s'', ''\\'', filesep);\n', local_matlab_string(manifest_file));
fprintf(fid, 'config_mat = strrep(''%s'', ''\\'', filesep);\n', local_matlab_string(config_mat));
fprintf(fid, 'load(config_mat, ''config'');\n');
fprintf(fid, '%% (Redundant addpath in case canlab_root differs between generation and run.)\n');
fprintf(fid, 'if isfield(config, ''canlab_root'') && ~isempty(config.canlab_root)\n');
fprintf(fid, '    addpath(genpath(config.canlab_root));\n');
fprintf(fid, 'end\n');
fprintf(fid, 'if isfield(config, ''extra_matlab_paths'')\n');
fprintf(fid, '    for path_idx = 1:numel(config.extra_matlab_paths)\n');
fprintf(fid, '        extra_path = local_cell_at(config.extra_matlab_paths, path_idx);\n');
fprintf(fid, '        if ~isempty(extra_path), addpath(genpath(strrep(extra_path, ''\\'', filesep))); end\n');
fprintf(fid, '    end\n');
fprintf(fid, 'end\n');
fprintf(fid, 'task_id = str2double(getenv(''SLURM_ARRAY_TASK_ID''));\n');
fprintf(fid, 'if isnan(task_id) || task_id < 1, task_id = str2double(getenv(''TASK_ID'')); end\n');
fprintf(fid, 'if isnan(task_id) || task_id < 1, task_id = 1; end\n');
fprintf(fid, 'if isfield(config, ''fmri_files'') && isfield(config, ''output_prefixes'')\n');
fprintf(fid, '    n_tasks = numel(config.fmri_files);\n');
fprintf(fid, '    if task_id > n_tasks, error(''Task id %%d exceeds configured tasks %%d.'', task_id, n_tasks); end\n');
fprintf(fid, '    subject = local_cell_at(config.subject_ids, task_id);\n');
fprintf(fid, '    run_label = local_cell_at(config.run_labels, task_id);\n');
fprintf(fid, '    fmri_file = strrep(local_cell_at(config.fmri_files, task_id), ''\\'', filesep);\n');
fprintf(fid, '    events_file = strrep(local_cell_at(config.events_files, task_id), ''\\'', filesep);\n');
fprintf(fid, '    output_prefix = strrep(local_cell_at(config.output_prefixes, task_id), ''\\'', filesep);\n');
fprintf(fid, '    output_mat = strrep(local_cell_at(config.output_mats, task_id), ''\\'', filesep);\n');
fprintf(fid, '    if isfield(config, ''spm_files'') && numel(config.spm_files) >= task_id, spm_file = strrep(local_cell_at(config.spm_files, task_id), ''\\'', filesep); else, spm_file = ''''; end\n');
fprintf(fid, '    if isfield(config, ''spm_runs'') && numel(config.spm_runs) >= task_id, spm_run = config.spm_runs(task_id); else, spm_run = 1; end\n');
fprintf(fid, 'else\n');
fprintf(fid, '    [subject, run_label, fmri_file, events_file, output_prefix, output_mat, spm_file, spm_run, n_tasks] = local_manifest_row(manifest_file, task_id);\n');
fprintf(fid, 'end\n');
fprintf(fid, 'if isempty(run_label), run_label = sprintf(''task-%%03d'', task_id); end\n');
fprintf(fid, 'fprintf(''Running HRF task %%d/%%d: %%s | %%s\\n'', task_id, n_tasks, subject, run_label);\n');
fprintf(fid, 'args = [config.pipeline_args, {''WriteWholeBrain'', true, ''WholeBrainOutputPrefix'', output_prefix, ''OutputMat'', output_mat}];\n');
fprintf(fid, 'if exist(''spm_file'', ''var'') && ~isempty(spm_file)\n');
fprintf(fid, '    args = [args, {''WholeBrainSPM'', spm_file, ''WholeBrainSPMRun'', spm_run}];\n');
fprintf(fid, '    fprintf(''  Tier B GKWY from SPM.mat: %%s (run %%d)\\n'', spm_file, spm_run);\n');
fprintf(fid, 'end\n');
fprintf(fid, 'results = run_hrf_pipeline(fmri_file, events_file, args{:});\n');
fprintf(fid, 'has_atlas_objs = isfield(config, ''atlas_objs'') && ~isempty(config.atlas_objs);\n');
fprintf(fid, 'if ~isempty(config.signature_sets) || ~isempty(config.image_sets) || ~isempty(config.atlas_obj) || has_atlas_objs\n');
fprintf(fid, '    wholebrain_models = local_wholebrain_models(results);\n');
fprintf(fid, '    for mi = 1:numel(wholebrain_models)\n');
fprintf(fid, '        model_name = wholebrain_models(mi).name;\n');
fprintf(fid, '        model_stats = wholebrain_models(mi).stats;\n');
fprintf(fid, '        score_context = sprintf(''task=%%d; subject=%%s; run=%%s'', task_id, subject, run_label);\n');
fprintf(fid, '        try\n');
fprintf(fid, '            score_status = hrf_score_one_prefix(output_prefix, ...\n');
fprintf(fid, '                ''ModelName'', model_name, ...\n');
fprintf(fid, '                ''NumModels'', numel(wholebrain_models), ...\n');
fprintf(fid, '                ''ScoreObjects'', config.score_objects, ...\n');
fprintf(fid, '                ''SignatureSets'', config.signature_sets, ...\n');
fprintf(fid, '                ''ImageSets'', config.image_sets, ...\n');
fprintf(fid, '                ''AtlasObj'', config.atlas_obj, ...\n');
fprintf(fid, '                ''AtlasName'', config.atlas_name, ...\n');
fprintf(fid, '                ''Regions'', config.regions, ...\n');
fprintf(fid, '                ''AtlasObjs'', local_config_field(config, ''atlas_objs'', {}), ...\n');
fprintf(fid, '                ''AtlasNames'', local_config_field(config, ''atlas_names'', {}), ...\n');
fprintf(fid, '                ''AtlasRegions'', local_config_field(config, ''atlas_regions'', {}), ...\n');
fprintf(fid, '                ''Normalize'', config.normalize, ...\n');
fprintf(fid, '                ''SimilarityMetric'', config.similarity_metric, ...\n');
fprintf(fid, '                ''StatsInput'', model_stats, ...\n');
fprintf(fid, '                ''Overwrite'', true, ...\n');
fprintf(fid, '                ''WarningContext'', score_context);\n');
fprintf(fid, '            for ei = 1:numel(score_status.errors)\n');
fprintf(fid, '                warning(''hrf_slurm_worker:ScoreError'', ''Score failure for model=%%s, object=%%s: %%s'', model_name, score_status.errors(ei).object, score_status.errors(ei).message);\n');
fprintf(fid, '            end\n');
fprintf(fid, '        catch score_err\n');
fprintf(fid, '            warning(''hrf_slurm_worker:ScoreFailed'', ''Scoring failed for model=%%s (task %%d): %%s'', model_name, task_id, score_err.message);\n');
fprintf(fid, '        end\n');
fprintf(fid, '    end\n');
fprintf(fid, 'end\n');
fprintf(fid, 'if config.make_animations\n');
fprintf(fid, '    movie_file = sprintf(''%%s_%%s_montage.mp4'', output_prefix, config.animation_object);\n');
fprintf(fid, '    wholebrain_models = local_wholebrain_models(results);\n');
fprintf(fid, '    hrf_animate_wholebrain_stats(wholebrain_models(1).stats, ''Object'', config.animation_object, ...\n');
fprintf(fid, '        ''Condition'', config.animation_condition, ''FrameRate'', config.animation_frame_rate, ...\n');
fprintf(fid, '        ''OutputFile'', movie_file);\n');
fprintf(fid, 'end\n');
fprintf(fid, '\n');
fprintf(fid, 'function models = local_wholebrain_models(results)\n');
fprintf(fid, 'models = struct(''name'', {}, ''stats'', {});\n');
fprintf(fid, 'if isfield(results, ''wholebrain_by_model'') && ~isempty(results.wholebrain_by_model)\n');
fprintf(fid, '    fields = fieldnames(results.wholebrain_by_model);\n');
fprintf(fid, '    for ii = 1:numel(fields)\n');
fprintf(fid, '        models(end + 1) = struct(''name'', fields{ii}, ''stats'', results.wholebrain_by_model.(fields{ii})); %%#ok<AGROW>\n');
fprintf(fid, '    end\n');
fprintf(fid, 'elseif isfield(results, ''wholebrain'')\n');
fprintf(fid, '    models = struct(''name'', ''wholebrain'', ''stats'', results.wholebrain);\n');
fprintf(fid, 'else\n');
fprintf(fid, '    error(''run_hrf_pipeline did not return wholebrain outputs.'');\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function val = local_cell_at(values, idx)\n');
fprintf(fid, 'if iscell(values), val = values{idx}; else, val = values(idx); end\n');
fprintf(fid, 'val = local_cell_to_char(val);\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function val = local_config_field(config, name, default_val)\n');
fprintf(fid, 'if isfield(config, name), val = config.(name); else, val = default_val; end\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function [subject, run_label, fmri_file, events_file, output_prefix, output_mat, spm_file, spm_run, n_tasks] = local_manifest_row(manifest_file, task_id)\n');
fprintf(fid, 'raw = local_read_manifest(manifest_file);\n');
fprintf(fid, 'if isempty(raw) || size(raw, 2) < 5, error(''Manifest must have at least 5 columns.''); end\n');
fprintf(fid, 'data_start = 1 + local_has_manifest_header(raw(1, :));\n');
fprintf(fid, 'n_tasks = size(raw, 1) - data_start + 1;\n');
fprintf(fid, 'if task_id > n_tasks, error(''Task id %%d exceeds manifest rows %%d.'', task_id, n_tasks); end\n');
fprintf(fid, 'row_idx = data_start + task_id - 1;\n');
fprintf(fid, 'headers = local_manifest_headers(raw, data_start);\n');
fprintf(fid, 'subject = local_manifest_value(raw, headers, row_idx, ''subject'', 2);\n');
fprintf(fid, 'run_label = local_manifest_value(raw, headers, row_idx, ''run_label'', 3);\n');
fprintf(fid, 'fmri_file = strrep(local_manifest_value(raw, headers, row_idx, ''fmri_file'', 4), ''\\'', filesep);\n');
fprintf(fid, 'events_file = strrep(local_manifest_value(raw, headers, row_idx, ''events_file'', 5), ''\\'', filesep);\n');
fprintf(fid, 'output_prefix = strrep(local_manifest_value(raw, headers, row_idx, ''output_prefix'', 6), ''\\'', filesep);\n');
fprintf(fid, 'output_mat = strrep(local_manifest_value(raw, headers, row_idx, ''output_mat'', 7), ''\\'', filesep);\n');
fprintf(fid, 'spm_file = strrep(local_manifest_value(raw, headers, row_idx, ''spm_file'', 8), ''\\'', filesep);\n');
fprintf(fid, 'spm_run = str2double(local_manifest_value(raw, headers, row_idx, ''spm_run'', 9));\n');
fprintf(fid, 'if isnan(spm_run) || spm_run < 1, spm_run = 1; end\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function headers = local_manifest_headers(raw, data_start)\n');
fprintf(fid, 'headers = containers.Map(''KeyType'', ''char'', ''ValueType'', ''double'');\n');
fprintf(fid, 'if data_start == 1, return; end\n');
fprintf(fid, 'for ii = 1:size(raw, 2)\n');
fprintf(fid, '    key = lower(strtrim(local_cell_to_char(raw{1, ii})));\n');
fprintf(fid, '    if ~isempty(key), headers(key) = ii; end\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function val = local_manifest_value(raw, headers, row_idx, name, fallback_col)\n');
fprintf(fid, 'if isKey(headers, name), col = headers(name); else, col = fallback_col; end\n');
fprintf(fid, 'if col > size(raw, 2), val = ''''; else, val = local_cell_to_char(raw{row_idx, col}); end\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function raw = local_read_manifest(manifest_file)\n');
fprintf(fid, 'try\n');
fprintf(fid, '    raw = readcell(manifest_file, ''Delimiter'', '','');\n');
fprintf(fid, 'catch\n');
fprintf(fid, '    fid2 = fopen(manifest_file, ''r'');\n');
fprintf(fid, '    if fid2 == -1, error(''Could not open manifest: %%s'', manifest_file); end\n');
fprintf(fid, '    c2 = onCleanup(@() fclose(fid2));\n');
fprintf(fid, '    lines = textscan(fid2, ''%%s'', ''Delimiter'', ''\\n'', ''Whitespace'', '''');\n');
fprintf(fid, '    lines = lines{1};\n');
fprintf(fid, '    raw = cell(numel(lines), 0);\n');
fprintf(fid, '    for ii = 1:numel(lines)\n');
fprintf(fid, '        parts = strsplit(lines{ii}, '','');\n');
fprintf(fid, '        raw(ii, 1:numel(parts)) = parts;\n');
fprintf(fid, '    end\n');
fprintf(fid, 'end\n');
fprintf(fid, 'raw = local_drop_empty_rows(raw);\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function raw = local_drop_empty_rows(raw)\n');
fprintf(fid, 'keep = true(size(raw, 1), 1);\n');
fprintf(fid, 'for ii = 1:size(raw, 1)\n');
fprintf(fid, '    row_text = strings(1, size(raw, 2));\n');
fprintf(fid, '    for jj = 1:size(raw, 2), row_text(jj) = string(local_cell_to_char(raw{ii, jj})); end\n');
fprintf(fid, '    keep(ii) = any(strlength(strtrim(row_text)) > 0);\n');
fprintf(fid, 'end\n');
fprintf(fid, 'raw = raw(keep, :);\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function tf = local_has_manifest_header(row)\n');
fprintf(fid, 'tokens = strings(1, numel(row));\n');
fprintf(fid, 'for ii = 1:numel(row), tokens(ii) = lower(strtrim(string(local_cell_to_char(row{ii})))); end\n');
fprintf(fid, 'tf = any(tokens == "subject") || any(tokens == "fmri_file") || any(tokens == "events_file") || any(tokens == "output_prefix");\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function s = local_cell_to_char(x)\n');
fprintf(fid, 'if iscell(x), x = x{1}; end\n');
fprintf(fid, 'try\n');
fprintf(fid, '    if ismissing(x), s = ''''; return; end\n');
fprintf(fid, 'catch\n');
fprintf(fid, 'end\n');
fprintf(fid, 'if isstring(x), s = char(x); elseif ischar(x), s = x; elseif isnumeric(x), s = num2str(x); else, s = char(string(x)); end\n');
fprintf(fid, 'end\n');
end

function local_write_sbatch(script_file, worker_file, n_tasks, log_dir, opts)
fid = fopen(script_file, 'w');
if fid == -1, error('Could not write SLURM script: %s', script_file); end
c = onCleanup(@() fclose(fid));

fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '#SBATCH --job-name=%s\n', char(opts.JobName));
fprintf(fid, '#SBATCH --array=1-%d\n', n_tasks);
fprintf(fid, '#SBATCH --time=%s\n', char(opts.Time));
fprintf(fid, '#SBATCH --mem=%s\n', char(opts.Mem));
fprintf(fid, '#SBATCH --cpus-per-task=%d\n', opts.CpusPerTask);
fprintf(fid, '#SBATCH --output=%s\n', local_posix_join(log_dir, '%x_%A_%a.out'));
fprintf(fid, '#SBATCH --error=%s\n', local_posix_join(log_dir, '%x_%A_%a.err'));

if ~isempty(char(opts.Partition))
    fprintf(fid, '#SBATCH --partition=%s\n', char(opts.Partition));
end
if ~isempty(char(opts.Account))
    fprintf(fid, '#SBATCH --account=%s\n', char(opts.Account));
end

extra = local_to_cell(opts.ExtraSlurmDirectives);
for i = 1:numel(extra)
    line = char(extra{i});
    if startsWith(strtrim(line), '#SBATCH')
        fprintf(fid, '%s\n', line);
    else
        fprintf(fid, '#SBATCH %s\n', line);
    end
end

fprintf(fid, '# Set SLURM ARRAY TASK ID if debugging interactively\n');
fprintf(fid, '\nset -euo pipefail\n\n');
fprintf(fid, 'export SLURM_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"\n\n');
if ~isempty(char(opts.ModuleLoad))
    fprintf(fid, 'module load %s\n\n', char(opts.ModuleLoad));
end
fprintf(fid, '%s -nodisplay -nosplash -batch "run(''%s'')"\n', ...
    char(opts.MatlabCommand), local_bash_matlab_path(worker_file));
end

function root = local_default_canlab_root()
% This file lives at CanlabCore/HRF_Est_Toolbox4/HRF_Toolbox_Pipeline/.
% Walk up three levels to reach the CanlabCore root (where @fmri_data etc live).
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

function local_ensure_parent(path_in)
parent_dir = fileparts(path_in);
if ~isempty(parent_dir) && ~exist(parent_dir, 'dir')
    mkdir(parent_dir);
end
end

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

function label = local_file_label(subject_id)
label = regexprep(char(subject_id), '[^\w.-]', '_');
end

function out = local_matlab_string(path_in)
out = strrep(char(path_in), '''', '''''');
end

function out = local_bash_matlab_path(path_in)
out = strrep(local_matlab_string(path_in), '\', '/');
end

function out = local_posix_join(folder, file_name)
out = strrep(fullfile(folder, file_name), '\', '/');
end
