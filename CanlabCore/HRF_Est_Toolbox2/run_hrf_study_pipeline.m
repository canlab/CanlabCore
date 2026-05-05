function study = run_hrf_study_pipeline(fmri_files, events_files, subject_ids, varargin)
%RUN_HRF_STUDY_PIPELINE Run run_hrf_pipeline across subjects/images.

if ischar(fmri_files) || isstring(fmri_files), fmri_files = cellstr(fmri_files); end
if ischar(events_files) || isstring(events_files), events_files = cellstr(events_files); end
if nargin < 3 || isempty(subject_ids)
    subject_ids = arrayfun(@(i) sprintf('sub-%03d', i), 1:numel(fmri_files), 'UniformOutput', false);
end
if isstring(subject_ids), subject_ids = cellstr(subject_ids); end

n = numel(fmri_files);
results = cell(1, n);
for i = 1:n
    results{i} = run_hrf_pipeline(fmri_files{i}, events_files{i}, varargin{:});
end

study = struct();
study.subject_ids = subject_ids;
study.results = results;
study.summary = hrf_summarize_study(results, subject_ids);

% Subject x time matrix of mean extracted timeseries
study.mean_timeseries = cell2mat(cellfun(@(r) r.timeseries(:)' , results, 'UniformOutput', false)');
end
