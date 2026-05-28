function obj = cat(varargin)
% cat - Concatenate fmri_hrf objects across (subject, run) axes.
%
% Usage
% -----
%   Hstudy = cat(1, Hb_sub01_run01, Hb_sub01_run02, ..., Hb_subN_runM)
%
% Concatenates the underlying fmri_data voxel-by-volume by delegating to
% fmri_data/cat. Then aligns HRF metadata: requires all inputs to share
% model_name, TR, and condition set (so the resulting stacked metadata
% table is interpretable). Stacks per-object (subject, run_label) into the
% metadata_table to preserve the source axis.
%
% v0 ignores the dimension argument and always concatenates along the
% volume/image axis. Multi-dim concatenation will come later.

if nargin < 2
    error('fmri_hrf:cat:NotEnoughArgs', 'cat requires at least two arguments.');
end

if isnumeric(varargin{1}) && isscalar(varargin{1})
    dim = varargin{1};
    objs = varargin(2:end);
    if dim ~= 1
        warning('fmri_hrf:cat:DimIgnored', ...
            'fmri_hrf/cat v0 only supports dim=1 (volume axis). Got dim=%d; treating as 1.', dim);
    end
else
    objs = varargin;
end

objs = objs(~cellfun(@isempty, objs));
if isempty(objs)
    obj = fmri_hrf();
    return
end

local_assert_alignable(objs);

% Delegate voxel-axis concatenation to fmri_data/cat.
parent_objs = cellfun(@(x) fmri_data(x), objs, 'UniformOutput', false);
parent_cat = cat(1, parent_objs{:});

% Build the stacked HRF metadata table by appending each input's
% metadata_table and tagging with subject + run_label.
meta_chunks = cell(numel(objs), 1);
for i = 1:numel(objs)
    M = objs{i}.metadata_table;
    if isempty(M)
        continue
    end
    M.subject = repmat(string(objs{i}.subject), height(M), 1);
    M.run_label = repmat(string(objs{i}.run_label), height(M), 1);
    meta_chunks{i} = M;
end
meta_chunks = meta_chunks(~cellfun(@isempty, meta_chunks));
if isempty(meta_chunks)
    stacked_meta = table();
else
    stacked_meta = vertcat(meta_chunks{:});
end

% Lift the concatenated fmri_data into a new fmri_hrf.
obj = fmri_hrf(parent_cat, ...
    'MetadataTable', stacked_meta, ...
    'ModelName', objs{1}.model_name, ...
    'Conditions', objs{1}.conditions, ...
    'TR', objs{1}.TR);
end


function local_assert_alignable(objs)
% All objects must share model_name, TR, conditions for concatenation to
% be meaningful. Subject/run can differ (that's the point).
ref = objs{1};
for i = 2:numel(objs)
    if ~strcmp(objs{i}.model_name, ref.model_name)
        error('fmri_hrf:cat:ModelMismatch', ...
            'Cannot concatenate fmri_hrf with different model_name: %s vs %s.', ...
            ref.model_name, objs{i}.model_name);
    end
    if ~isequaln(objs{i}.TR, ref.TR)
        error('fmri_hrf:cat:TRMismatch', ...
            'Cannot concatenate fmri_hrf with different TR: %g vs %g.', ...
            ref.TR, objs{i}.TR);
    end
    if ~isequal(sort(cellstr(string(objs{i}.conditions))), sort(cellstr(string(ref.conditions))))
        error('fmri_hrf:cat:ConditionMismatch', ...
            'Cannot concatenate fmri_hrf with different condition sets.');
    end
end
end
