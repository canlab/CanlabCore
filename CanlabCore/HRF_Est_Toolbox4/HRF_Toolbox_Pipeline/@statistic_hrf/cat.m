function obj = cat(varargin)
% cat - Concatenate statistic_hrf objects across (subject, run) axes.
%
% Usage
% -----
%   Tstudy = cat(1, Ht_sub01_run01, Ht_sub01_run02, ..., Ht_subN_runM)
%
% Mirrors fmri_hrf/cat: delegates voxel-axis concatenation to
% statistic_image/cat and stacks the HRF metadata tables, tagging each
% chunk with its subject and run_label.

if nargin < 2
    error('statistic_hrf:cat:NotEnoughArgs', 'cat requires at least two arguments.');
end

if isnumeric(varargin{1}) && isscalar(varargin{1})
    dim = varargin{1};
    objs = varargin(2:end);
    if dim ~= 1
        warning('statistic_hrf:cat:DimIgnored', ...
            'statistic_hrf/cat v0 only supports dim=1. Got dim=%d; treating as 1.', dim);
    end
else
    objs = varargin;
end

objs = objs(~cellfun(@isempty, objs));
if isempty(objs)
    obj = statistic_hrf();
    return
end

local_assert_alignable(objs);

parent_objs = cellfun(@(x) statistic_image(x), objs, 'UniformOutput', false);
parent_cat = cat(1, parent_objs{:});

meta_chunks = cell(numel(objs), 1);
for i = 1:numel(objs)
    M = objs{i}.metadata_table;
    if isempty(M), continue; end
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

obj = statistic_hrf(parent_cat, ...
    'MetadataTable', stacked_meta, ...
    'ModelName', objs{1}.model_name, ...
    'Conditions', objs{1}.conditions, ...
    'TR', objs{1}.TR);
end


function local_assert_alignable(objs)
ref = objs{1};
for i = 2:numel(objs)
    if ~strcmp(objs{i}.model_name, ref.model_name)
        error('statistic_hrf:cat:ModelMismatch', ...
            'Cannot concatenate statistic_hrf with different model_name: %s vs %s.', ...
            ref.model_name, objs{i}.model_name);
    end
    if ~isequaln(objs{i}.TR, ref.TR)
        error('statistic_hrf:cat:TRMismatch', ...
            'Cannot concatenate statistic_hrf with different TR: %g vs %g.', ...
            ref.TR, objs{i}.TR);
    end
    if ~isequal(sort(cellstr(string(objs{i}.conditions))), sort(cellstr(string(ref.conditions))))
        error('statistic_hrf:cat:ConditionMismatch', ...
            'Cannot concatenate statistic_hrf with different condition sets.');
    end
end
end
