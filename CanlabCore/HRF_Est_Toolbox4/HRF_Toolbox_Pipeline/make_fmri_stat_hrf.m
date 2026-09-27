function [Hb, Ht] = make_fmri_stat_hrf(source, varargin)
% make_fmri_stat_hrf - Build a paired (fmri_hrf, statistic_hrf) from a single fit.
%
% This is the canonical entry point to the new HRF object classes. The two
% sibling classes (fmri_hrf for beta/amplitude, statistic_hrf for t/p/SE)
% are designed to be constructed together so their HRF metadata is
% guaranteed to be aligned.
%
% Usage
% -----
%   1) From an in-memory wholebrain_by_model struct:
%      [Hb, Ht] = make_fmri_stat_hrf(results.wholebrain_by_model.sfir, ...
%                                    'Subject', 'sub-01', ...
%                                    'RunLabel', 'task-pain_run-01', ...
%                                    'ModelName', 'sfir', ...
%                                    'TR', 0.8)
%
%   2) From a NIfTI prefix (loads beta + t + se + metadata from disk):
%      [Hb, Ht] = make_fmri_stat_hrf('/path/to/sub-01_hrf_sfir', ...
%                                    'Subject', 'sub-01', ...
%                                    'ModelName', 'sfir', ...
%                                    'TR', 0.8)
%
%   3) From a result.mat path:
%      [Hb, Ht] = make_fmri_stat_hrf('/path/to/sub-01_hrf_results.mat', ...
%                                    'ModelName', 'sfir', ...
%                                    'Subject', 'sub-01')
%
% Inputs
% ------
%   source     The first arg dispatches on type:
%               struct  -> in-memory wholebrain struct with .b, .t, fields
%                          (and optionally .ste / .se for SE)
%               char/string ending in .mat -> result.mat path
%               char/string (any other)    -> NIfTI prefix; reads
%                          <prefix>_beta.nii, _t.nii, _se.nii, _metadata.csv
%
% Required name-value pairs
% -------------------------
%   'Subject', 'ModelName', 'TR'.  (RunLabel optional but recommended.)
%
% Optional name-value pairs
% -------------------------
%   'MetadataTable'   - override the metadata table resolution.
%   'Conditions'      - condition labels; derived from metadata if not given.
%   'DesignMatrix'    - design matrix (needed for residuals later).
%   'DesignInfo'      - struct describing design columns.
%   'NoVerbose'       - suppress fmri_data load chatter (default true).
%
% Returns
% -------
%   Hb   fmri_hrf, beta side.
%   Ht   statistic_hrf, t side. Empty if no t image is resolvable from the source.
%
% See also: fmri_hrf, statistic_hrf, hrf_fit_wholebrain_stats.

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('source');
p.addParameter('Subject', '', @(x) ischar(x) || isstring(x));
p.addParameter('RunLabel', '', @(x) ischar(x) || isstring(x));
p.addParameter('ModelName', '', @(x) ischar(x) || isstring(x));
p.addParameter('TR', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Conditions', {}, @(x) iscell(x) || isstring(x));
p.addParameter('MetadataTable', table(), @(x) isempty(x) || istable(x));
p.addParameter('DesignMatrix', [], @(x) isnumeric(x) || isempty(x));
p.addParameter('DesignInfo', struct(), @isstruct);
p.addParameter('NoVerbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(source, varargin{:});
opts = p.Results;

[beta_obj, t_obj, se_obj, metadata, source_paths] = local_resolve_source(source, opts);

if isempty(opts.MetadataTable) && ~isempty(metadata)
    opts.MetadataTable = metadata;
end

nv_pairs = local_pack_nv(opts, source_paths);

% Build the beta side.
Hb = fmri_hrf(beta_obj, nv_pairs{:});

% Build the t side if a t image is available; else derive from beta + se if both.
if ~isempty(t_obj)
    Ht = statistic_hrf(t_obj, nv_pairs{:});
elseif ~isempty(se_obj)
    Ht = to_statistic_hrf(Hb, 'SE', se_obj);
else
    Ht = statistic_hrf();
end
end


% =========================================================================
% Source dispatch
% =========================================================================
function [beta_obj, t_obj, se_obj, metadata, source_paths] = local_resolve_source(source, opts)
beta_obj = [];
t_obj = [];
se_obj = [];
metadata = table();
source_paths = struct();

if isstruct(source)
    [beta_obj, t_obj, se_obj, metadata, source_paths] = local_from_wholebrain_struct(source);
    return
end

if ischar(source) || isstring(source)
    source = char(source);
    if endsWith(lower(source), '.mat')
        [beta_obj, t_obj, se_obj, metadata, source_paths] = local_from_result_mat(source, opts);
    else
        [beta_obj, t_obj, se_obj, metadata, source_paths] = local_from_prefix(source, opts);
    end
    return
end

error('make_fmri_stat_hrf:UnknownSource', ...
    'source must be a wholebrain struct, a result.mat path, or a NIfTI prefix. Got: %s.', ...
    class(source));
end

function [beta_obj, t_obj, se_obj, metadata, source_paths] = local_from_wholebrain_struct(s)
beta_obj = [];
t_obj = [];
se_obj = [];
metadata = table();
source_paths = struct();

if isfield(s, 'b'), beta_obj = s.b; end
if isfield(s, 'beta'), beta_obj = s.beta; end
if isfield(s, 't'), t_obj = s.t; end
if isfield(s, 'tstat'), t_obj = s.tstat; end
if isfield(s, 'se'), se_obj = s.se; end
if isfield(s, 'ste'), se_obj = s.ste; end

if isempty(se_obj) && ~isempty(beta_obj) && isprop(beta_obj, 'ste') && ~isempty(beta_obj.ste)
    se_obj = beta_obj;
    se_obj.dat = beta_obj.ste;
end

if isfield(s, 'metadata_table') && istable(s.metadata_table)
    metadata = s.metadata_table;
elseif isfield(s, 'metadata') && istable(s.metadata)
    metadata = s.metadata;
end
end

function [beta_obj, t_obj, se_obj, metadata, source_paths] = local_from_prefix(prefix, opts)
load_args = {};
if logical(opts.NoVerbose), load_args = {'noverbose'}; end

beta_file = [prefix '_beta.nii'];
t_file = [prefix '_t.nii'];
se_file = [prefix '_se.nii'];
meta_file = [prefix '_metadata.csv'];

source_paths = struct('beta', beta_file, 't', t_file, 'se', se_file, 'metadata', meta_file);

beta_obj = [];
if exist(beta_file, 'file') == 2
    beta_obj = fmri_data(beta_file, load_args{:});
end

t_obj = [];
if exist(t_file, 'file') == 2
    t_obj = statistic_image(fmri_data(t_file, load_args{:}));
    t_obj.type = 'T';
end

se_obj = [];
if exist(se_file, 'file') == 2
    se_obj = statistic_image(fmri_data(se_file, load_args{:}));
    se_obj.type = 'HRF beta standard error';
end

metadata = table();
if exist(meta_file, 'file') == 2
    try
        metadata = readtable(meta_file, 'TextType', 'string');
    catch
    end
end

if isempty(beta_obj)
    error('make_fmri_stat_hrf:MissingBeta', ...
        'Could not find %s. Check the prefix.', beta_file);
end
end

function [beta_obj, t_obj, se_obj, metadata, source_paths] = local_from_result_mat(mat_file, opts)
beta_obj = [];
t_obj = [];
se_obj = [];
metadata = table();
source_paths = struct('result_mat', mat_file);

if exist(mat_file, 'file') ~= 2
    error('make_fmri_stat_hrf:MissingResultMat', 'result.mat not found: %s', mat_file);
end

S = load(mat_file, 'results');
if ~isfield(S, 'results')
    error('make_fmri_stat_hrf:NoResultsField', ...
        'result.mat does not contain a top-level "results" struct.');
end

model_name = char(opts.ModelName);
if isempty(model_name)
    error('make_fmri_stat_hrf:NeedModelName', ...
        'Provide ''ModelName'' when constructing from a result.mat.');
end

R = S.results;
model_field = matlab.lang.makeValidName(lower(model_name));
if ~isfield(R, 'wholebrain_by_model') || ~isfield(R.wholebrain_by_model, model_field)
    error('make_fmri_stat_hrf:MissingModelInResultMat', ...
        'result.mat has no wholebrain_by_model.%s entry.', model_field);
end

[beta_obj, t_obj, se_obj, metadata, ~] = local_from_wholebrain_struct(R.wholebrain_by_model.(model_field));

if isempty(metadata) && isfield(R, 'wholebrain_metadata_by_model') && isfield(R.wholebrain_metadata_by_model, model_field)
    metadata = R.wholebrain_metadata_by_model.(model_field);
end
end

function nv = local_pack_nv(opts, source_paths)
nv = {};
add = @(name, val) [nv, {name, val}]; %#ok<NASGU>
nv = [nv, {'Subject', char(opts.Subject)}];
nv = [nv, {'RunLabel', char(opts.RunLabel)}];
nv = [nv, {'ModelName', char(opts.ModelName)}];
nv = [nv, {'TR', opts.TR}];
if ~isempty(opts.Conditions), nv = [nv, {'Conditions', opts.Conditions}]; end
if ~isempty(opts.MetadataTable), nv = [nv, {'MetadataTable', opts.MetadataTable}]; end
if ~isempty(opts.DesignMatrix), nv = [nv, {'DesignMatrix', opts.DesignMatrix}]; end
if ~isempty(fieldnames(opts.DesignInfo)), nv = [nv, {'DesignInfo', opts.DesignInfo}]; end
if ~isempty(fieldnames(source_paths)), nv = [nv, {'SourcePaths', source_paths}]; end
end
