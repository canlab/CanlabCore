% statistic_hrf: HRF inferential container (t / p / SE), paired with fmri_hrf.
%
% Subclass of statistic_image. Holds whole-brain HRF t-maps (or p, SE) from
% a single fit with the metadata to interpret them: condition labels, lag
% indices, TR. All threshold/multi_threshold/orthviews/convert2mask methods
% are inherited from statistic_image and work unchanged.
%
% Multi-model studies use one statistic_hrf per (subject, run, model),
% paired with the matching fmri_hrf. See make_fmri_stat_hrf for the
% constructor that builds them in lockstep.
%
% Usage
% -----
%   Ht = statistic_hrf()                                       % empty
%   Ht = statistic_hrf(statistic_image_obj, 'Name', Val, ...)  % lift
%
% Required name-value pairs (same as fmri_hrf): MetadataTable, Subject,
% RunLabel, ModelName, TR. Optional: Conditions, DesignMatrix, DesignInfo,
% SourcePaths.
%
% See also: fmri_hrf, make_fmri_stat_hrf, statistic_image, threshold.

classdef statistic_hrf < statistic_image
    properties
        hrf_meta_version = 1
        metadata_table = table()
        subject = ''
        run_label = ''
        model_name = ''
        conditions = {}
        TR = NaN
        design_matrix = []
        design_info = struct()
        source_paths = struct()
    end

    methods
        function obj = statistic_hrf(varargin)
            lifting = nargin > 0 && (isa(varargin{1}, 'statistic_image') || isa(varargin{1}, 'image_vector') || isa(varargin{1}, 'fmri_data'));
            if lifting
                parent_args = {};
                nv_args = varargin(2:end);
            else
                parent_args = varargin;
                nv_args = {};
            end

            obj = obj@statistic_image(parent_args{:});

            if lifting
                src = varargin{1};
                obj = local_copy_parent_fields(obj, src);
            end

            if ~isempty(nv_args)
                obj = local_apply_hrf_metadata(obj, nv_args);
            end

            if nargin > 0
                local_validate(obj);
            end
        end

        function disp(obj)
            if isempty(obj.dat)
                fprintf('  statistic_hrf  (empty)\n');
                return
            end
            n_vox = size(obj.dat, 1);
            n_vol = size(obj.dat, 2);
            n_cond = numel(obj.conditions);
            n_lag = NaN;
            if ~isempty(obj.metadata_table) && any(strcmp('lag_index', obj.metadata_table.Properties.VariableNames))
                n_lag = numel(unique(obj.metadata_table.lag_index));
            end
            fprintf('  statistic_hrf  subject=%s  run=%s  model=%s  type=%s\n', ...
                local_disp(obj.subject), local_disp(obj.run_label), local_disp(obj.model_name), local_disp(obj.type));
            fprintf('                 %d voxels x %d volumes  (%d conditions x %d lags)  TR=%.3gs\n', ...
                n_vox, n_vol, n_cond, n_lag, obj.TR);
            if ~isempty(obj.conditions)
                fprintf('                 conditions: %s\n', strjoin(cellstr(string(obj.conditions)), ', '));
            end
        end
    end
end


% =========================================================================
% Local helpers
% =========================================================================
function obj = local_copy_parent_fields(obj, src)
hrf_only = {'hrf_meta_version', 'metadata_table', 'subject', 'run_label', ...
    'model_name', 'conditions', 'TR', 'design_matrix', 'design_info', ...
    'source_paths'};

src_props = properties(src);
for i = 1:numel(src_props)
    f = src_props{i};
    if ismember(f, hrf_only), continue; end
    if isprop(obj, f)
        try
            obj.(f) = src.(f);
        catch
        end
    end
end
end

function obj = local_apply_hrf_metadata(obj, nv_args)
p = inputParser;
p.KeepUnmatched = false;
p.addParameter('MetadataTable', table(), @(x) isempty(x) || istable(x));
p.addParameter('Subject', '', @(x) ischar(x) || isstring(x));
p.addParameter('RunLabel', '', @(x) ischar(x) || isstring(x));
p.addParameter('ModelName', '', @(x) ischar(x) || isstring(x));
p.addParameter('Conditions', {}, @(x) iscell(x) || isstring(x));
p.addParameter('TR', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('DesignMatrix', [], @(x) isnumeric(x) || isempty(x));
p.addParameter('DesignInfo', struct(), @isstruct);
p.addParameter('SourcePaths', struct(), @isstruct);
p.parse(nv_args{:});
opts = p.Results;

if ~isempty(opts.MetadataTable), obj.metadata_table = opts.MetadataTable; end
if ~isempty(char(opts.Subject)), obj.subject = char(opts.Subject); end
if ~isempty(char(opts.RunLabel)), obj.run_label = char(opts.RunLabel); end
if ~isempty(char(opts.ModelName)), obj.model_name = char(opts.ModelName); end
if ~isempty(opts.Conditions)
    obj.conditions = cellstr(string(opts.Conditions));
elseif ~isempty(obj.metadata_table) && any(strcmp('condition', obj.metadata_table.Properties.VariableNames))
    obj.conditions = unique(cellstr(string(obj.metadata_table.condition)), 'stable');
end
if ~isnan(opts.TR), obj.TR = opts.TR; end
if ~isempty(opts.DesignMatrix), obj.design_matrix = opts.DesignMatrix; end
if ~isempty(fieldnames(opts.DesignInfo)), obj.design_info = opts.DesignInfo; end
if ~isempty(fieldnames(opts.SourcePaths)), obj.source_paths = opts.SourcePaths; end
end

function local_validate(obj)
if isempty(obj.dat), return; end
n_vol = size(obj.dat, 2);
if ~isempty(obj.metadata_table) && height(obj.metadata_table) ~= n_vol
    error('statistic_hrf:MetadataVolumeMismatch', ...
        'metadata_table has %d rows but underlying statistic_image has %d volumes.', ...
        height(obj.metadata_table), n_vol);
end
end

function s = local_disp(x)
if isempty(x)
    s = '<unset>';
else
    s = char(x);
end
end
