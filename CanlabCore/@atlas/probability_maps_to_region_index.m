function dat = probability_maps_to_region_index(dat)
% probability_maps_to_region_index Rebuild integer index labels from probability_maps.
%
% Use dat.probability_maps to rebuild the integer vector of index labels
% in dat.dat: each voxel is assigned to the region with maximum
% probability, with voxels having all-zero or all-NaN probabilities set
% to 0. If some regions are missing (no nonzero column in
% probability_maps), they and their corresponding labels and
% label_descriptions are dropped from the atlas, and the remaining index
% values are renumbered consecutively.
%
% :Usage:
% ::
%
%     dat = probability_maps_to_region_index(dat)
%
% :Inputs:
%
%   **dat:**
%        An atlas-class object with a populated probability_maps field.
%
% :Outputs:
%
%   **dat:**
%        Atlas-class object whose .dat is an integer index vector
%        consistent with the columns of probability_maps.
%
% :Notes:
%
% This script should not invoke remove_empty, replace_empty, or
% resampling because it is invoked on atlas objects that may have
% mismatched probability_maps and dat fields, causing such functions to
% fail.
%
% :See also:
%   - check_properties
%   - num_regions
%   - select_atlas_subset

% Start: dat has one image per region, with probability values
% convert to integer vector

[maxval, condf] = max(double(full(dat.probability_maps)),[], 2);   % double is safer

allempty = maxval == 0 | isnan(maxval);

condf(allempty) = 0;

dat.dat = int32(condf);
n_regions = num_regions(dat);
n_maps = find(any(full(dat.probability_maps) > 0,1));
if length(n_maps) < n_regions
    dropped_ind = find(~ismember(1:n_regions,n_maps));
    
    % note, labels must be last in this list or else it would end up being
    % asigned to the labels variable which is needed subsequently
    fnames = {'probability_maps', 'labels_2','labels_3','labels_4',...
        'labels_5','label_descriptions','labels'};
    assert(strcmp(fnames{end},'labels'));
    for f = 1:length(fnames)            
        if length(dat.labels) ~= n_regions
            for i = 1:length(dropped_ind)
                warning('Dropping %s index %d', fnames{f}, dropped_ind(i));
            end
        end
        if length(dat.(fnames{f})) == n_regions
            labels = dat.(fnames{f});
            dat.(fnames{f})(dropped_ind) = [];
        elseif size(dat.(fnames{f}),2) == n_regions
            dat.(fnames{f})(:,dropped_ind) = [];
        end
    end
    if length(dat.labels) == n_regions
        for i = 1:length(dropped_ind)
            warning('Dropping region %d: %s', dropped_ind(i), labels{dropped_ind(i)});
        end
    end
    if ismember(0,unique(dat.dat))
        [~,~, dat.dat] = unique(dat.dat);
        dat.dat = dat.dat - 1;
    else
        [~,~, dat.dat] = unique(dat.dat);
    end
end


