function study = run_hrf_study_pipeline(fmri_files, events_files, subject_ids, varargin)
%RUN_HRF_STUDY_PIPELINE Run run_hrf_pipeline across subjects/images.

if ischar(fmri_files) || isstring(fmri_files), fmri_files = cellstr(fmri_files); end
if ischar(events_files) || isstring(events_files), events_files = cellstr(events_files); end
if nargin < 3 || isempty(subject_ids)
    subject_ids = arrayfun(@(i) sprintf('sub-%03d', i), 1:numel(fmri_files), 'UniformOutput', false);
end
if isstring(subject_ids), subject_ids = cellstr(subject_ids); end

[study_opts, pipeline_args] = local_parse_study_options(varargin{:});

n = numel(fmri_files);
if numel(events_files) ~= n
    error('fmri_files and events_files must contain the same number of files.');
end
if numel(subject_ids) ~= n
    error('subject_ids must be empty or contain one id per fMRI file.');
end
results = cell(1, n);
errors = cell(1, n);

if ~isempty(study_opts.wholebrain_output_dir) && ~exist(study_opts.wholebrain_output_dir, 'dir')
    mkdir(study_opts.wholebrain_output_dir);
end

wholebrain_output_dir = study_opts.wholebrain_output_dir;

if study_opts.use_parallel_subjects
    if isempty(gcp('nocreate'))
        parpool;
    end
    parfor i = 1:n
        [results{i}, errors{i}] = local_run_one(fmri_files{i}, events_files{i}, subject_ids{i}, pipeline_args, wholebrain_output_dir);
    end
else
    for i = 1:n
        [results{i}, errors{i}] = local_run_one(fmri_files{i}, events_files{i}, subject_ids{i}, pipeline_args, wholebrain_output_dir);
    end
end

failed = ~cellfun(@isempty, errors);
if any(failed)
    for i = find(failed)
        warning('Problem with %s: %s', fmri_files{i}, errors{i});
    end
    if ~study_opts.continue_on_error
        error('run_hrf_study_pipeline failed for %d file(s). First error: %s', sum(failed), errors{find(failed, 1)});
    end
end

success = ~cellfun(@isempty, results);
valid_results = results(success);
valid_subject_ids = subject_ids(success);

study = struct();
study.subject_ids = subject_ids;
study.results = results;
study.success = success;
study.errors = errors;
study.summary = hrf_summarize_study(valid_results, valid_subject_ids);

% Subject x time matrix of mean extracted timeseries
if isempty(valid_results)
    study.mean_timeseries = [];
else
    ts = cellfun(@(r) r.timeseries(:)', valid_results, 'UniformOutput', false);
    ts_len = cellfun(@numel, ts);
    if all(ts_len == ts_len(1))
        study.mean_timeseries = cell2mat(ts');
    else
        study.mean_timeseries = ts;
    end
end
end

function [opts, pipeline_args] = local_parse_study_options(varargin)
opts = struct('use_parallel_subjects', false, 'continue_on_error', true, 'wholebrain_output_dir', '');
pipeline_args = varargin;
i = 1;
while i <= numel(pipeline_args)
    if ischar(pipeline_args{i}) || isstring(pipeline_args{i})
        key = lower(char(pipeline_args{i}));
        switch key
            case 'useparallelsubjects'
                opts.use_parallel_subjects = logical(pipeline_args{i + 1});
                pipeline_args(i:i+1) = [];
                continue
            case 'continueonerror'
                opts.continue_on_error = logical(pipeline_args{i + 1});
                pipeline_args(i:i+1) = [];
                continue
            case 'wholebrainoutputdir'
                opts.wholebrain_output_dir = char(pipeline_args{i + 1});
                pipeline_args(i:i+1) = [];
                continue
        end
    end
    i = i + 1;
end
end

function [result, err_msg] = local_run_one(fmri_file, events_file, subject_id, pipeline_args, wholebrain_output_dir)
result = [];
err_msg = '';
try
    args = pipeline_args;
    if ~isempty(wholebrain_output_dir)
        prefix = fullfile(wholebrain_output_dir, [local_file_label(subject_id) '_hrf']);
        args = [args, {'WholeBrainOutputPrefix', prefix}];
    end
    result = run_hrf_pipeline(fmri_file, events_file, args{:});
catch ME
    err_msg = ME.message;
end
end

function label = local_file_label(subject_id)
label = regexprep(char(subject_id), '[^\w.-]', '_');
end
