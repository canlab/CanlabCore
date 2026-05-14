function T = hrf_collect_wholebrain_outputs(root_dir, varargin)
%HRF_COLLECT_WHOLEBRAIN_OUTPUTS Index HRF 4D outputs for second-level work.
%
% T = hrf_collect_wholebrain_outputs(root_dir)
%
% Looks recursively for files written by hrf_fit_wholebrain_stats:
%   *_beta.nii
%   *_t.nii
%   *_se.nii
%   *_p.nii
%   *_t_thresh.nii
%   *_metadata.csv
% and optional map-score files written by hrf_apply_maps_to_wholebrain:
%   *_beta_map_scores.csv
%   *_t_map_scores.csv
% and optional result MAT files written by run_hrf_pipeline:
%   *_results.mat

p = inputParser;
p.addRequired('root_dir', @(x) ischar(x) || isstring(x));
p.addParameter('OutputCsv', '', @(x) ischar(x) || isstring(x));
p.parse(root_dir, varargin{:});
opts = p.Results;

root_dir = char(root_dir);
beta_files = dir(fullfile(root_dir, '**', '*_beta.nii'));

rows = {};
for i = 1:numel(beta_files)
    beta_path = fullfile(beta_files(i).folder, beta_files(i).name);
    prefix = erase(beta_path, '_beta.nii');
    t_path = [prefix '_t.nii'];
    se_path = [prefix '_se.nii'];
    p_path = [prefix '_p.nii'];
    t_thresh_path = [prefix '_t_thresh.nii'];
    metadata_path = [prefix '_metadata.csv'];
    beta_scores_path = [prefix '_beta_map_scores.csv'];
    t_scores_path = [prefix '_t_map_scores.csv'];
    result_mat_path = [prefix '_results.mat'];

    [~, prefix_name] = fileparts(prefix);
    subject = local_subject_from_name(prefix_name);

    rows(end + 1, :) = {subject, prefix, beta_path, local_existing(t_path), ...
        local_existing(se_path), local_existing(p_path), ...
        local_existing(t_thresh_path), local_existing(metadata_path), ...
        local_existing(beta_scores_path), local_existing(t_scores_path), ...
        local_existing(result_mat_path)}; %#ok<AGROW>
end

var_names = {'subject', 'prefix', 'beta_file', 't_file', 'se_file', 'p_file', 'thresholded_t_file', ...
    'metadata_file', 'beta_scores_file', 't_scores_file', 'result_mat_file'};
if isempty(rows)
    T = cell2table(cell(0, numel(var_names)), 'VariableNames', var_names);
else
    T = cell2table(rows, 'VariableNames', var_names);
end

if ~isempty(opts.OutputCsv)
    writetable(T, char(opts.OutputCsv));
end
end

function out = local_existing(path_in)
if exist(path_in, 'file')
    out = path_in;
else
    out = '';
end
end

function subject = local_subject_from_name(name)
tok = regexp(name, '(sub-[A-Za-z0-9]+)', 'tokens', 'once');
if isempty(tok)
    subject = name;
else
    subject = tok{1};
end
end
