function avg = hrf_pooled_wholebrain_animation(dirs, varargin)
%HRF_POOLED_WHOLEBRAIN_ANIMATION Pool whole-brain HRF maps across dirs and animate.
%
% Convenience wrapper that pools the whole-brain HRF maps of one OR MORE study
% output directories and writes group (and optionally per-subject) montage
% animations over the peristimulus lags, for the conditions you specify.
%
% It collects each directory with hrf_collect_wholebrain_outputs, concatenates
% the resulting tables (the SAME subject id appearing in two dirs has that
% subject's runs pooled -- e.g. both bodysite dirs of distractmap), and hands
% the combined table to hrf_make_average_montage_animations, which does the
% per-condition run->subject->group averaging and renders the movies.
%
% NOTE: the 'Conditions' you pass here are the FIT conditions stored in the
% map metadata (e.g. 'rest_stim_ttl_1', 'nback-stimblock_ttl_1' for
% distractmap), NOT the raw BIDS event trial_types used for the Granger
% condition segmentation. Glob wildcards are allowed.
%
% :Usage:
% ::
%     dirs = {'...\hrf_outputs_lf_distractmap', '...\hrf_outputs_obs_distractmap'};
%     avg = hrf_pooled_wholebrain_animation(dirs, ...
%         'Model','sfir', 'Object','beta', ...
%         'Conditions', {'rest_stim_ttl_1','nback-stimblock_ttl_1'}, ...
%         'OutputDir', '...\pooled_distractmap_anim', ...
%         'GroupStatistic','t', 'GroupCorrection','fdr', 'GroupAlpha',0.05, ...
%         'MakeSubjectAnimations', false);
%
% :Inputs:
%   **dirs:** a study output directory (char) OR a cell array of directories
%             whose maps are pooled. Each must contain the whole-brain map
%             NIfTIs + *_metadata.csv that hrf_collect_wholebrain_outputs reads.
%
% :Optional Inputs:
%   All name-value arguments are passed through to
%   hrf_make_average_montage_animations. The ones you will most often set:
%     'Model'        - 'sfir' | 'canonical' | 'spline' | 'fir' (default '' =
%                      all models present -- usually set this).
%     'Object'       - 'beta' (default) | 't'.
%     'Conditions'   - cellstr of fit-condition labels/globs to animate.
%                      Default {} = every condition found (one movie each).
%     'OutputDir'    - where to write the .mp4 files (required to get movies).
%     'GroupStatistic' - 'mean' (default) or 't' (one-sample across subjects).
%     'GroupCorrection'/'GroupAlpha' - 'none'|'fdr', 0.05 (for the 't' movie).
%     'MakeSubjectAnimations' / 'MakeGroupAnimations' - default true/true.
%     'FPS','FrameStep','ColorMap','ColorLimits','Threshold' - rendering.
%   **'Verbose':** print the pooling summary (default true).
%
% :Output:
%   **avg:** the struct returned by hrf_make_average_montage_animations, plus
%            .pooled_dirs and .input_table (the concatenated collection table).
%
% See also: hrf_make_average_montage_animations, hrf_collect_wholebrain_outputs,
%           hrf_make_montage_animation, hrf_causality.

p = inputParser;
p.KeepUnmatched = true;       % everything else flows to the averager
p.addRequired('dirs', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(dirs, varargin{:});
verbose = logical(p.Results.Verbose);

dl = local_dir_list(dirs);
IT = table();
for i = 1:numel(dl)
    if exist(fullfile(dl{i}), 'dir') ~= 7
        warning('hrf_pooled_wholebrain_animation:NoDir', 'Directory not found, skipping: %s', dl{i});
        continue
    end
    Ti = hrf_collect_wholebrain_outputs(dl{i});
    if isempty(Ti) || height(Ti) == 0
        warning('hrf_pooled_wholebrain_animation:NoMaps', 'No whole-brain maps collected in: %s', dl{i});
        continue
    end
    Ti.source_dir = repmat(string(dl{i}), height(Ti), 1);
    IT = local_vcat(IT, Ti);
end
if isempty(IT) || height(IT) == 0
    error('hrf_pooled_wholebrain_animation:Empty', 'No whole-brain records pooled from the %d dir(s).', numel(dl));
end

if verbose
    nsubj = numel(unique(cellstr(string(IT.subject))));
    fprintf('hrf_pooled_wholebrain_animation: pooled %d dir(s) -> %d records, %d subjects\n', ...
        numel(dl), height(IT), nsubj);
end

% Forward the remaining name-value pairs (minus our own 'Verbose') untouched.
fwd = local_strip_param(varargin, 'Verbose');
avg = hrf_make_average_montage_animations(IT, fwd{:});
avg.pooled_dirs = dl(:)';
avg.input_table = IT;
end


% =========================================================================
function L = local_dir_list(x)
if iscell(x)
    L = cellfun(@(c) char(string(c)), x, 'uni', 0);
elseif isstring(x) && ~isscalar(x)
    L = cellstr(x);
else
    L = {char(string(x))};
end
end


function T = local_vcat(T, Ti)
% Concatenate two collection tables on their COMMON columns (defensive against
% a dir that has a slightly different column set).
if isempty(T) || height(T) == 0
    T = Ti;
    return
end
common = intersect(T.Properties.VariableNames, Ti.Properties.VariableNames, 'stable');
T = [T(:, common); Ti(:, common)];
end


function args = local_strip_param(args, name)
% Remove a name-value pair (case-insensitive) from a varargin cell list.
keep = true(1, numel(args));
i = 1;
while i <= numel(args) - 1
    if (ischar(args{i}) || isstring(args{i})) && strcmpi(char(args{i}), name)
        keep(i) = false; keep(i + 1) = false;
        i = i + 2;
    else
        i = i + 1;
    end
end
args = args(keep);
end
