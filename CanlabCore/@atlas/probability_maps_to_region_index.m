function dat = probability_maps_to_region_index(dat)
% Use dat.probability_maps to rebuild integer vector of index labels (dat.dat)
%
% dat = probability_maps_to_region_index(dat)
%
% Note: this script should not invoke remove_empty, replace_empty or
% resampling because it is invoked on atlas objects that may have
% mismatched probability_maps and dat fields, causing such functions to 
% fail.

% Start: dat has one image per region, with probability values
% convert to integer vector

[maxval, condf] = max(double(full(dat.probability_maps)),[], 2);   % double is safer

allempty = maxval == 0 | isnan(maxval);

condf(allempty) = 0;

dat.dat = int32(condf);
n_regions = num_regions(dat);
% I chagned this line sunday night, and also the n_maps reference below
n_maps = find(any(full(dat.probability_maps) > 0,1));
if length(n_maps) < n_regions
    dropped_ind = find(~ismember(1:n_regions,n_maps));
    
    % note, labels must be last in this list or lese it would end up being
    % asigned to the labels variable which is needed subsequently
    fnames = {'probability_maps', 'labels_2','labels_3','labels_4',...
        'labels_5','label_descriptions','labels'};
    assert(strcmp(fnames{end},'labels'));
    for f = 1:length(fnames)            
        if length(dat.labels) ~= n_regions
            for i = 1:length(dropped_ind)
                warning('Dropping %s index %d', fnames{f}, dropped_ind);
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
    if ismember(dat.dat, 0)
        [~,~, dat.dat] = unique(dat.dat);
        dat.dat = dat.dat - 1;
    else
        [~,~, dat.dat] = unique(dat.dat);
    end
end


