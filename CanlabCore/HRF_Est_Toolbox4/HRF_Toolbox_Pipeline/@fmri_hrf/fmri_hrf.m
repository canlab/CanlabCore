% fmri_hrf: HRF amplitude/beta container, paired with statistic_hrf for inference.
%
% Subclass of fmri_data. Holds whole-brain HRF beta maps from a single fit
% (one subject, one run, one model) along with the metadata needed to
% interpret them: condition labels, lag indices, TR, design matrix.
%
% Multi-model studies use one fmri_hrf per (subject, run, model). Cross-model
% comparison (Phase 4) is handled by a separate hrf_modelset container, not
% by stuffing multiple models into one fmri_hrf object.
%
% Usage
% -----
%   Hb = fmri_hrf()                                    % empty constructor
%   Hb = fmri_hrf(fmri_data_obj, 'Name', Val, ...)     % lift an fmri_data
%
% The multi-source dispatch (build from a wholebrain_by_model struct, a
% result.mat path, or a NIfTI prefix) lives in make_fmri_stat_hrf, which
% calls this constructor with an already-loaded fmri_data plus metadata.
%
% Required name-value pairs when constructing from an fmri_data:
%   'MetadataTable'  - table with one row per 4D volume. Must include columns
%                      condition, lag_index, lag_seconds, image_label.
%   'Subject'        - char/string, BIDS subject id.
%   'RunLabel'       - char/string, BIDS task-run id.
%   'ModelName'      - 'fir' | 'sfir' | 'canonical' | 'spline'.
%   'TR'             - scalar seconds.
%
% Optional name-value pairs:
%   'Conditions'     - cellstr of condition labels (derived from metadata if empty).
%   'DesignMatrix'   - design matrix (volumes x regressors), for residuals().
%   'DesignInfo'     - struct with column labels / FIR layout.
%   'SourcePaths'    - struct of file paths the maps came from.
%   'Timeseries'     - optional struct array per signature/ROI. See note below.
%
% Inherited from fmri_data
% ------------------------
% All voxel storage, masking, reading/writing, resample, etc. are inherited
% unchanged. Methods overridden on fmri_hrf:
%   cat       - aligns HRF metadata across (subject, run) before delegating
%               to fmri_data/cat for the underlying voxel concatenation.
%
% Notes
% -----
% .timeseries is an optional struct array storing 1D signal extractions
% (signature time series, ROI averages) that share the subject/run/condition
% metadata of the voxel data but have completely different storage (time x
% signatures rather than voxels x volumes). v0 only stores it; methods that
% operate on the .timeseries side will be added when needed.
%
% See also: statistic_hrf, make_fmri_stat_hrf, fmri_data, hrf_apply_maps_to_wholebrain.

classdef fmri_hrf < fmri_data
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
        timeseries = []
    end

    methods
        function obj = fmri_hrf(varargin)
            % Always call the parent constructor. If we are lifting an
            % existing fmri_data, call parent with no args and copy fields
            % after; otherwise pass through.
            lifting = nargin > 0 && (isa(varargin{1}, 'fmri_data') || isa(varargin{1}, 'image_vector'));
            if lifting
                parent_args = {};
                nv_args = varargin(2:end);
            else
                parent_args = varargin;
                nv_args = {};
            end

            obj = obj@fmri_data(parent_args{:});

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
            % HRF-aware display summary.
            if isempty(obj.dat)
                fprintf('  fmri_hrf  (empty)\n');
                return
            end
            n_vox = size(obj.dat, 1);
            n_vol = size(obj.dat, 2);
            n_cond = numel(obj.conditions);
            n_lag = NaN;
            if ~isempty(obj.metadata_table) && any(strcmp('lag_index', obj.metadata_table.Properties.VariableNames))
                n_lag = numel(unique(obj.metadata_table.lag_index));
            end
            fprintf('  fmri_hrf  subject=%s  run=%s  model=%s\n', ...
                local_disp(obj.subject), local_disp(obj.run_label), local_disp(obj.model_name));
            fprintf('            %d voxels x %d volumes  (%d conditions x %d lags)  TR=%.3gs\n', ...
                n_vox, n_vol, n_cond, n_lag, obj.TR);
            if ~isempty(obj.conditions)
                fprintf('            conditions: %s\n', strjoin(cellstr(string(obj.conditions)), ', '));
            end
            if ~isempty(obj.timeseries)
                fprintf('            timeseries: %d entries attached\n', numel(obj.timeseries));
            end
        end
    end
end


% =========================================================================
% Local helpers
% =========================================================================
function obj = local_copy_parent_fields(obj, src)
% Copy data + metadata fields from src (fmri_data or image_vector) into obj.
% Skips properties that fmri_hrf overrides.
hrf_only = {'hrf_meta_version', 'metadata_table', 'subject', 'run_label', ...
    'model_name', 'conditions', 'TR', 'design_matrix', 'design_info', ...
    'source_paths', 'timeseries'};

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
p.addParameter('Timeseries', [], @(x) isempty(x) || isstruct(x));
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
if ~isempty(opts.Timeseries), obj.timeseries = opts.Timeseries; end
end

function local_validate(obj)
if isempty(obj.dat), return; end
n_vol = size(obj.dat, 2);
if ~isempty(obj.metadata_table) && height(obj.metadata_table) ~= n_vol
    error('fmri_hrf:MetadataVolumeMismatch', ...
        'metadata_table has %d rows but underlying fmri_data has %d volumes.', ...
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
