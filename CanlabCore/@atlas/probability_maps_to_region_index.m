function dat = probability_maps_to_region_index(dat)
% Use dat.probability_maps to rebuild integer vector of index labels (dat.dat)
%
% dat = probability_maps_to_region_index(dat)
%

% Start: dat has one image per region, with probability values
% convert to integer vector

[maxval, condf] = max(double(full(dat.probability_maps)),[], 2);   % double is safer

allempty = all(dat.probability_maps == 0, 2) | isnan(maxval);  % some out-of-mask values may get NaNs

condf(allempty) = 0;

dat.dat = int32(condf);
n_regions = size(dat.probability_maps,2);
if length(unique(condf)) < n_regions
    dropped_ind = find(~ismember(1:n_regions,unique(condf)));
    labels = cell(length(dropped_ind),1);
    fnames = {'labels_2','labels_3','labels_4','labels_5','label_descriptions','labels'};
    for f = 1:length(fnames)
        if length(dat.(fnames{f})) == n_regions
            labels = dat.(fnames{f});
            dat.(fnames{f})(dropped_ind) = [];
        end
    end
    for i = 1:length(dropped_ind)
        warning('Dropping region %d: %s', dropped_ind(i), labels{dropped_ind(i)});
    end
end


