function obj = cat(varargin)
% cat - Concatenate fmri_hrf objects across (subject, run) axes.
%
% Usage
% -----
%   Hstudy = cat(1, Hb_sub01_run01, Hb_sub01_run02, ..., Hb_subN_runM)
%
% Stacks the underlying voxel .dat matrices along the volume axis and
% stacks each input's HRF metadata_table, tagging chunks with subject and
% run_label so the source axis is preserved. The output's other parent
% fields (volInfo, mask, etc.) are copied from the first input.
%
% v0 ignores the dimension argument and always concatenates along the
% volume axis. Voxel count must match across inputs.

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

% Stack the .dat matrices along the volume axis.
base = objs{1};
n_vox = size(base.dat, 1);
all_dat = base.dat;
for i = 2:numel(objs)
    if size(objs{i}.dat, 1) ~= n_vox
        error('fmri_hrf:cat:VoxelMismatch', ...
            'Cannot cat: input %d has %d voxels, expected %d.', ...
            i, size(objs{i}.dat, 1), n_vox);
    end
    all_dat = [all_dat, objs{i}.dat]; %#ok<AGROW>
end

% Build the stacked HRF metadata table, tagging chunks with subject + run_label.
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

% Build the result by value-copying the first input and updating data + meta.
obj = base;
obj.dat = all_dat;
obj.metadata_table = stacked_meta;
% Reset image-axis flags to match new volume count.
n_vol = size(all_dat, 2);
if isprop(obj, 'removed_images')
    obj.removed_images = false(n_vol, 1);
end
end


function local_assert_alignable(objs)
% All objects must share model_name, TR, conditions for concatenation to
% be meaningful. Subject and run_label can differ -- that's the point.
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
