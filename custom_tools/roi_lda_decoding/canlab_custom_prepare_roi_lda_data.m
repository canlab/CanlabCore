function prep = canlab_custom_prepare_roi_lda_data(C_obj, Group_labels, atlas_input, varargin)
% Prepare subject-by-voxel ROI data for ROI-wise decoding.
%
% prep = canlab_custom_prepare_roi_lda_data(C_obj, Group_labels, atlas_input)
%
% Inputs
% -------------------------------------------------------------------------
% C_obj         1 x nStudies cell array of contrast-like Canlab objects.
%               Each cell must contain a .dat matrix with either:
%               - subjects x voxels, or
%               - voxels x subjects
%
% Group_labels  nStudies labels that assign each study to a decoding class.
%
% atlas_input   Atlas object or atlas keyword accepted by load_atlas.
%
% Optional name-value pairs
% -------------------------------------------------------------------------
% 'verbose'     true/false, default true
%
% Output
% -------------------------------------------------------------------------
% prep is a struct with aligned subject-by-voxel data and ROI metadata.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'C_obj', @(x) iscell(x) && ~isempty(x));
addRequired(p, 'Group_labels', @(x) numel(x) == numel(C_obj));
addRequired(p, 'atlas_input', @(x) ischar(x) || isstring(x) || isa(x, 'atlas'));
addParameter(p, 'verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, C_obj, Group_labels, atlas_input, varargin{:});

verbose = p.Results.verbose;
nStudies = numel(C_obj);

if ischar(atlas_input) || isstring(atlas_input)
    atlas_obj = load_atlas(char(atlas_input));
else
    atlas_obj = atlas_input;
end

atlas_obj = replace_empty(atlas_obj);
ref_obj = replace_empty(C_obj{1});

isdiff = compare_space(ref_obj, atlas_obj);
if isdiff == 2
    error('%s: atlas or reference object is missing volInfo.', mfilename);
elseif ~(isdiff == 0 || isdiff == 3)
    if verbose
        fprintf('%s: resampling atlas to study space.\n', mfilename);
    end
    atlas_obj = resample_space(atlas_obj, ref_obj);
end

study_data = cell(nStudies, 1);
study_nsubjects = zeros(nStudies, 1);
shared_mask = logical(atlas_obj.volInfo.image_indx);

for i = 1:nStudies
    obj = replace_empty(C_obj{i});
    isdiff = compare_space(obj, ref_obj);

    if isdiff == 2
        error('%s: C_obj{%d} is missing volInfo.', mfilename, i);
    elseif ~(isdiff == 0 || isdiff == 3)
        error('%s: all C_obj entries must already be in the same voxel space.', mfilename);
    end

    study_data{i} = local_subject_by_voxel(obj);
    study_nsubjects(i) = size(study_data{i}, 1);
    shared_mask = shared_mask & logical(obj.volInfo.image_indx);
end

atlas_vec = local_inmask_vector(atlas_obj);
atlas_keep = shared_mask(atlas_obj.volInfo.wh_inmask);
atlas_index = double(round(atlas_vec(atlas_keep)));

region_ids = unique(atlas_index);
region_ids(~isfinite(region_ids) | region_ids == 0) = [];

for i = 1:nStudies
    obj = replace_empty(C_obj{i});
    keep_mask = shared_mask(obj.volInfo.wh_inmask);
    study_data{i} = study_data{i}(:, keep_mask);
end

X = cat(1, study_data{:});
[study_labels_numeric, class_names] = local_encode_labels(Group_labels);
subject_labels = repelem(study_labels_numeric(:), study_nsubjects(:));
study_id_subject = repelem((1:nStudies)', study_nsubjects(:));

region_labels = cell(numel(region_ids), 1);
for i = 1:numel(region_ids)
    if ~isempty(atlas_obj.labels) && region_ids(i) <= numel(atlas_obj.labels)
        region_labels{i} = atlas_obj.labels{region_ids(i)};
    else
        region_labels{i} = sprintf('ROI_%d', region_ids(i));
    end
end

prep = struct();
prep.X = X;
prep.study_data = study_data;
prep.study_nsubjects = study_nsubjects;
prep.subject_labels = subject_labels;
prep.study_labels_numeric = study_labels_numeric(:);
prep.class_names = class_names;
prep.study_id_subject = study_id_subject;
prep.shared_mask = shared_mask;
prep.atlas_obj = atlas_obj;
prep.atlas_index = atlas_index;
prep.region_ids = region_ids(:);
prep.region_labels = region_labels;
prep.nsubjects_total = size(X, 1);
prep.nvoxels_shared = size(X, 2);

end

function X = local_subject_by_voxel(obj)
nvox = numel(obj.volInfo.wh_inmask);

if size(obj.dat, 2) == nvox
    X = double(obj.dat);
elseif size(obj.dat, 1) == nvox
    X = double(obj.dat');
else
    error('%s: could not infer subject/voxel orientation for object of class %s.', ...
        mfilename, class(obj));
end
end

function vec = local_inmask_vector(obj)
if size(obj.dat, 1) == numel(obj.volInfo.wh_inmask)
    vec = obj.dat(:, 1);
elseif size(obj.dat, 1) == numel(obj.volInfo.image_indx)
    vec = obj.dat(logical(obj.volInfo.image_indx), 1);
else
    error('%s: atlas .dat has an unexpected size.', mfilename);
end
end

function [labels_numeric, class_names] = local_encode_labels(labels_in)
if iscategorical(labels_in)
    labels_cell = cellstr(labels_in(:));
elseif isstring(labels_in)
    labels_cell = cellstr(labels_in(:));
elseif iscellstr(labels_in)
    labels_cell = labels_in(:);
elseif iscell(labels_in)
    labels_cell = cellfun(@local_to_char, labels_in(:), 'UniformOutput', false);
elseif isnumeric(labels_in) || islogical(labels_in)
    [u, ~, labels_numeric] = unique(labels_in(:), 'stable');
    class_names = cellstr(string(u));
    return
else
    error('%s: unsupported Group_labels type.', mfilename);
end

[class_names, ~, labels_numeric] = unique(labels_cell, 'stable');
end

function out = local_to_char(in)
if isstring(in)
    out = char(in);
elseif ischar(in)
    out = in;
elseif isnumeric(in) || islogical(in)
    out = char(string(in));
else
    error('%s: unsupported label element type.', mfilename);
end
end
