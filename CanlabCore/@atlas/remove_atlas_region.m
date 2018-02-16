function atlas_obj = remove_atlas_region(atlas_obj, varargin)
% Removes region(s) from atlas, by names or index values
% See select_atlas_subset for input format for region names/values
%
% atlas_obj = remove_atlas_region(atlas_obj, [region names cell or index values vector])
% 


[atlas_obj, has_pmaps] = check_properties(atlas_obj);

% Get list to_remove in index values
% ---------------------------------------------------------------------

[~, to_remove] = select_atlas_subset(atlas_obj, varargin{:});

% Remove from data
% Then rebuild index values, or adjust mask values
% ---------------------------------------------------------------------

if has_pmaps
    
    atlas_obj.probability_maps(:, to_remove) = [];
    
    atlas_obj = probability_maps_to_region_index(atlas_obj); % rebuild dat
    
else
    
    to_subtract = cumsum(to_remove);  % subtract
    
    for i = 1:length(to_subtract)
        
        wh = atlas_obj.dat == i;
        
        atlas_obj.dat(wh) = atlas_obj.dat(wh) - to_subtract(i);
        
    end
    
end

% Remove from labels and descriptions
% ---------------------------------------------------------------------

n = num_regions(atlas_obj);

myfields = {'labels' 'label_descriptions' 'labels_2' 'labels_3' 'labels_4' 'labels_5'};

for i = 1:length(myfields)
    
    if ~isempty(atlas_obj.(myfields{i})) && length(atlas_obj.(myfields{i})) == n
        
        atlas_obj.(myfields{i})(to_remove) = [];
        
    end
    
end

% these fields may be in char arrays, handle differently
myfields = {'references' 'image_names' 'fullpath'};

for i = 1:length(myfields)
    
    if ~isempty(atlas_obj.(myfields{i})) && size(atlas_obj.(myfields{i}), 1) == n
        
        atlas_obj.(myfields{i})(to_remove, :) = [];
        
    end
    
end


end % Function
