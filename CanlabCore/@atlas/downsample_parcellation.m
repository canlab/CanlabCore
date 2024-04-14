function new_atlas_obj = downsample_parcellation(obj, varargin)
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
%     varargin - either a vector of new labels to use, a string or leave empty. 
%                  If you provide a vector, its length should equal number of
%                  regions in original atlas, with redundant labels for regions
%                  you want to combine in the downsampled atlas.
%                  If you provide a string, it should be one of 
%                  {'labels_2', 'labels_3', 'labels_4', 'labels_5'} which will
%                  then be used as the downsampling vector.
%                  If you don't provide anything the contents of labels_2 is used
%                  by default.
%
% Output ::
%     atlas_obj - a new atlas object, reindexed according to labelfield

fprintf('Downsampling %s parcels\n', obj.atlas_name);

concat_desc = false;

if isempty(varargin)
    labelfield = 'labels_2';
else
    for i = 1:length(varargin)
        if ismember(varargin{i},{'labels_2','labels_3','labels_3','labels_4','labels_5'})
            labelfield=varargin{i};
        elseif strcmp(varargin{i},'concat_label_descriptions')
            concat_desc = true;
        elseif iscell(varargin{i})
            if length(varargin{i}) ~= num_regions(obj)
                error('New labels has length %d which does not match input atlas parcel count (%d)',length(varargin{1}), num_regions(obj));
            end
            labelfield='labels';
            obj.labels = varargin{1};
        else
            error('Unrecognized input datatype');
        end
    end
end

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
%kept_lbl_ind = find((n_unique_lbls <= n_unique_lbls(new_ind)) & ...
%    (n_lbls == length(obj.(labelfield))));

% find labels with valid fields
kept_lbl_ind = find((n_lbls == length(obj.(labelfield))));
% keep those with indices that are greater than or equal to the desired label
kept_lbl_ind = kept_lbl_ind(kept_lbl_ind >= new_ind);

rm_lbl_ind = find(~ismember(1:length(valid_lbls), 1:length(kept_lbl_ind)));

% merge high resolution regions
[new_lbl_parcels, lbl_exp, c] = unique(obj.(new_lbl),'stable');
new_parcels = cell(1,length(new_lbl_parcels));

% check that higher order labels are nested
for i = 1:length(lbl_exp) % number of unique labels
    this_ind = find(i == c);
    for j = 2:length(kept_lbl_ind)
        higher_order_lbls = obj.(valid_lbls{kept_lbl_ind(j)});
        uniq_ho_lbls = unique(higher_order_lbls(this_ind));
        if length(uniq_ho_lbls) > 1
            lbls_with_amp = strcat(uniq_ho_lbls,repmat({' & '},1,length(uniq_ho_lbls)));
            lbls_with_amp = cat(2,lbls_with_amp{:});
            lbls_with_amp = lbls_with_amp(1:end-3); % drop trailing ampersand
            obj.(valid_lbls{kept_lbl_ind(j)})(this_ind) = {lbls_with_amp};

            warning('Higher order labels are not nested. Merging select entries of %s into %s', valid_lbls{kept_lbl_ind(j)},lbls_with_amp);
        end
    end
end
        

% init progress watcher
n_reg = length(new_lbl_parcels);
last_msg = sprintf('Creating new region 1/%d',n_reg);
fprintf('%s',last_msg);

% working with sparse matrices is much faster for atlases with many
% parcels, while atlases with few parcels are likely to run quickly
% regardless
obj.probability_maps = sparse(double(obj.probability_maps));
for i = 1:length(new_lbl_parcels)
    % this can be parallelized but it's very memory intensive for some
    % atlases (e.g. canlab2023) and parallelization makes the problem worse
    % by several factors.

    % update progress watcher
    delete_last_chars = length(last_msg);
    fprintf('%s',char(8*ones(1,delete_last_chars)))
    new_msg = sprintf('Creating new region %d/%d',i,n_reg);
    fprintf('%s',new_msg);
    last_msg = new_msg;

    new_parcel{i} = obj.select_atlas_subset(new_lbl_parcels(i), new_lbl, 'exact', 'flatten');
    if concat_desc
        % append semicolon to all constituent descriptions
        label_desc = cellfun(@(x1)[x1, ';'],  obj.select_atlas_subset(new_lbl_parcels(i), new_lbl, 'exact').label_descriptions, 'UniformOutput', false);
        new_parcel{i}.label_descriptions = {strcat(label_desc{:})};
    end
end
fprintf('\n');

% init progress watcher
n_reg = length(new_parcel);
last_msg = sprintf('Merging new region 1/%d',n_reg);
fprintf('%s',last_msg);

% Combine parcels
new_atlas_obj = new_parcel{1};
for i = 2:length(new_parcel)
    % update progress watcher
    delete_last_chars = length(last_msg);
    fprintf('%s',char(8*ones(1,delete_last_chars)))
    new_msg = sprintf('Merging new region %d/%d',i,n_reg);
    fprintf('%s',new_msg);
    last_msg = new_msg;

    new_atlas_obj = new_atlas_obj.merge_atlases(new_parcel{i},'noverbose');
end
fprintf('\n');
new_atlas_obj.probability_maps = sparse(new_atlas_obj.probability_maps);

for i = 1:length(kept_lbl_ind)
    % increment by 1 
    old_lbls = obj.(valid_lbls{kept_lbl_ind(i)});

    new_atlas_obj.(valid_lbls{i}) = old_lbls(lbl_exp(1:length(lbl_exp)));
end

% delete finer parcellations
for ind = rm_lbl_ind(:)'
    new_atlas_obj.(valid_lbls{ind}) = {};
end

new_atlas_obj.references = unique(new_atlas_obj.references,'rows');

end