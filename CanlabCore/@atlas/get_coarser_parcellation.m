function new_atlas_obj = get_coarser_parcellation(obj, labelfield)
% An atlas_obj has multiple label fields meant to label nested parcellations
% at different levels of granularity. labels is the default which must
% contain unique entries and correspond to the number of unique indices
% in the (voxelwise) parcellation index map, but labels_2 ... labels_5
% may contain non-unique entries corresponding to coarser parcellations.
% This function remaps indices and probability maps to correspond to one
% of these coarser parcellations. It is a lossy operations and all labels
% of finer parcellation than the selected level will be discarded.
%
% Input ::
%     obj - an atlas_object. Must have at least one of labels_2 ...
%                 labels_5 defined
%     labelfield - which label field to substitute for labels. Must be one of
%                  {'labels_2', 'labels_3', 'labels_4', 'labels_5'}
%
% Output ::
%     atlas_obj - a new atlas object, reindexed according to labelfield

% input check
valid_lbls = {'labels', 'labels_2', 'labels_3', 'labels_4', 'labels_5'};
if ~ismember(labelfield, valid_lbls)
    error('labelfield must be one of {''labels_2'', ''labels_3'', ''labels_4'', ''labels_5''})');
end

% check that labels are ordered with finer labels first
n_unique_lbls = zeros(1,length(valid_lbls));
n_lbls = zeros(1,length(valid_lbls));
for i = 1:length(valid_lbls)
    lbl = valid_lbls{i};
    n_unique_lbls(i) = length(unique(obj.(lbl)));
    n_lbls(i) = length(obj.(lbl));
end

new_ind = find(strcmp(valid_lbls,labelfield));
if n_unique_lbls(new_ind) > n_unique_lbls(1)
    error('labelfield %s has %d labels, which is finer grained than labels which has %d unique labels',labelfield, n_unique_lbls(new_ind), n_unique_lbls(1));
end

new_lbl = valid_lbls{new_ind};
% note: we remove any label fields below if they're not of the appropriate
% length to correspond to the original parcellation
kept_lbl_ind = find((n_unique_lbls <= n_unique_lbls(new_ind)) & ...
    (n_lbls == length(obj.(labelfield))));

rm_lbl_ind = find(~ismember(1:length(valid_lbls), 1:length(kept_lbl_ind)));

[new_lbl_parcels, lbl_exp] = unique(obj.(new_lbl),'stable');
new_parcels = cell(1,length(new_lbl_parcels));
parfor i = 1:length(new_lbl_parcels)
    new_parcel{i} = obj.select_atlas_subset(new_lbl_parcels(i), new_lbl, 'exact', 'flatten');
end

new_atlas_obj = new_parcel{1};
for i = 2:length(new_parcel)
    new_atlas_obj = new_atlas_obj.merge_atlases(new_parcel{i});
end
for i = 1:length(kept_lbl_ind)
    % increment by 1 
    old_lbls = obj.(valid_lbls{kept_lbl_ind(i)});

    new_atlas_obj.(valid_lbls{i}) = old_lbls(lbl_exp(1:length(lbl_exp)));
end

% delete finer parcellations
for ind = rm_lbl_ind(:)'
    new_atlas_obj.(valid_lbls{ind}) = {};
end

end
