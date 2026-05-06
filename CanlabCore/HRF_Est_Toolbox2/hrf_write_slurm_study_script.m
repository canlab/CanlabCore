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
p.addParameter('CanlabRoot', local_default_canlab_root(), @(x) ischar(x) || isstring(x));
p.addParameter('PipelineArgs', {}, @(x) iscell(x));
p.addParameter('SignatureSets', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ImageSets', {}, @(x) ischar(x) || iscell(x) || isstring(x) || isa(x, 'image_vector'));
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

manifest = local_manifest_table(fmri_files, events_files, subject_ids, output_dir);
writetable(manifest, paths.manifest_file);

config = struct();
config.canlab_root = char(opts.CanlabRoot);
config.pipeline_args = opts.PipelineArgs;
config.signature_sets = local_to_cell(opts.SignatureSets);
config.image_sets = local_to_cell(opts.ImageSets);
config.score_objects = local_to_cell(opts.ScoreObjects);
config.similarity_metric = char(opts.SimilarityMetric);
config.make_animations = logical(opts.MakeAnimations);
config.animation_object = char(opts.AnimationObject);
config.animation_condition = char(opts.AnimationCondition);
config.animation_frame_rate = opts.AnimationFrameRate;
save(paths.config_mat, 'config');

local_write_worker(paths.worker_file, paths.manifest_file, paths.config_mat);
local_write_sbatch(paths.script_file, paths.worker_file, n, log_dir, opts);
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

function T = local_manifest_table(fmri_files, events_files, subject_ids, output_dir)
n = numel(fmri_files);
subject = cell(n, 1);
fmri_file = cell(n, 1);
events_file = cell(n, 1);
output_prefix = cell(n, 1);
output_mat = cell(n, 1);

for i = 1:n
    label = local_file_label(subject_ids{i});
    subject{i} = char(subject_ids{i});
    fmri_file{i} = char(fmri_files{i});
    events_file{i} = char(events_files{i});
    output_prefix{i} = fullfile(output_dir, [label '_hrf']);
    output_mat{i} = fullfile(output_dir, [label '_hrf_results.mat']);
end

T = table((1:n)', subject, fmri_file, events_file, output_prefix, output_mat, ...
    'VariableNames', {'index', 'subject', 'fmri_file', 'events_file', 'output_prefix', 'output_mat'});
end

function local_write_worker(worker_file, manifest_file, config_mat)
fid = fopen(worker_file, 'w');
if fid == -1, error('Could not write worker file: %s', worker_file); end
c = onCleanup(@() fclose(fid));

fprintf(fid, '%% Auto-generated by hrf_write_slurm_study_script.m\n');
fprintf(fid, 'manifest_file = strrep(''%s'', ''\\'', filesep);\n', local_matlab_string(manifest_file));
fprintf(fid, 'config_mat = strrep(''%s'', ''\\'', filesep);\n', local_matlab_string(config_mat));
fprintf(fid, 'load(config_mat, ''config'');\n');
fprintf(fid, 'addpath(genpath(config.canlab_root));\n');
fprintf(fid, 'manifest = readtable(manifest_file, ''TextType'', ''string'');\n');
fprintf(fid, 'task_id = str2double(getenv(''SLURM_ARRAY_TASK_ID''));\n');
fprintf(fid, 'if isnan(task_id) || task_id < 1, task_id = str2double(getenv(''TASK_ID'')); end\n');
fprintf(fid, 'if isnan(task_id) || task_id < 1, task_id = 1; end\n');
fprintf(fid, 'if task_id > height(manifest), error(''Task id %%d exceeds manifest rows %%d.'', task_id, height(manifest)); end\n');
fprintf(fid, 'row = manifest(task_id, :);\n');
fprintf(fid, 'fprintf(''Running HRF task %%d/%%d: %%s\\n'', task_id, height(manifest), char(row.subject));\n');
fprintf(fid, 'args = [config.pipeline_args, {''WriteWholeBrain'', true, ''WholeBrainOutputPrefix'', char(row.output_prefix), ''OutputMat'', char(row.output_mat)}];\n');
fprintf(fid, 'results = run_hrf_pipeline(char(row.fmri_file), char(row.events_file), args{:});\n');
fprintf(fid, 'if ~isempty(config.signature_sets) || ~isempty(config.image_sets)\n');
fprintf(fid, '    for oi = 1:numel(config.score_objects)\n');
fprintf(fid, '        object_name = char(config.score_objects{oi});\n');
fprintf(fid, '        score_csv = sprintf(''%%s_%%s_map_scores.csv'', char(row.output_prefix), object_name);\n');
fprintf(fid, '        hrf_apply_maps_to_wholebrain(results.wholebrain, ''Object'', object_name, ...\n');
fprintf(fid, '            ''SignatureSets'', config.signature_sets, ''ImageSets'', config.image_sets, ...\n');
fprintf(fid, '            ''SimilarityMetric'', config.similarity_metric, ''OutputCsv'', score_csv);\n');
fprintf(fid, '    end\n');
fprintf(fid, 'end\n');
fprintf(fid, 'if config.make_animations\n');
fprintf(fid, '    movie_file = sprintf(''%%s_%%s_montage.mp4'', char(row.output_prefix), config.animation_object);\n');
fprintf(fid, '    hrf_animate_wholebrain_stats(results.wholebrain, ''Object'', config.animation_object, ...\n');
fprintf(fid, '        ''Condition'', config.animation_condition, ''FrameRate'', config.animation_frame_rate, ...\n');
fprintf(fid, '        ''OutputFile'', movie_file);\n');
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

fprintf(fid, '\nset -euo pipefail\n\n');
if ~isempty(char(opts.ModuleLoad))
    fprintf(fid, 'module load %s\n\n', char(opts.ModuleLoad));
end
fprintf(fid, '%s -nodisplay -nosplash -batch "run(''%s'')"\n', ...
    char(opts.MatlabCommand), local_bash_matlab_path(worker_file));
end

function root = local_default_canlab_root()
root = fileparts(fileparts(mfilename('fullpath')));
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
