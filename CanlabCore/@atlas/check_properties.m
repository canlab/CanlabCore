function obj = check_properties(obj)
% Check properties and enforce some variable types
%
% obj = check_properties(obj)

n_regions = num_regions(obj);

% Data
has_pmaps = ~isempty(obj.probability_maps) && size(obj.probability_maps, 2) == n_regions;
has_index = ~isempty(obj.dat);

% Create index if needed
if has_pmaps && ~has_index
    obj = probability_maps_to_region_index(obj);
end

% could check for valid data here and fix, warn, or error

valid_pmaps = all(obj.probability_maps(:) >= 0) && all(obj.probability_maps(:) <= 1);
valid_index = size(obj.dat, 2) == 1 && all(abs(obj.dat(:) - round(obj.dat(:))) < 100 * eps);

n_regions_from_indx = length(unique(obj.dat(obj.dat ~= 0)));

if ~(max(obj.dat) == n_regions) || n_regions_from_indx ~= n_regions
    warning('Integer index in obj.dat should have one integer value per label, without skipping any integers!!');
end

% Enforce variable types
% ----------------------------------------------------------------------------

obj.probability_maps = sparse(double(obj.probability_maps));
obj.dat = int32(full(round(obj.dat)));

myfields = {'labels' 'labels_2' 'labels_3' 'labels_4' 'labels_5' 'property_descriptions', 'label_descriptions'};

for i = 1:length(myfields)
    
    if ~iscell(obj.(myfields{i})), obj.(myfields{i}) = {obj.(myfields{i})}; end
    
end

% Check for names/labels and create placeholders if missing
% ----------------------------------------------------------------------------

if isempty(obj.label_descriptions)
    obj.label_descriptions = cell(n_regions, 1);
end

if isempty(obj.labels)
    obj.labels = cell(1, n_regions);
end
    
for i = 1:n_regions
    
    if length(obj.labels) < i || isempty(obj.labels{i})
        
        obj.labels{i} = sprintf('Reg_%d', i);
        
    end
    
    if length(obj.label_descriptions) < i || isempty(obj.label_descriptions{i})
        
        obj.label_descriptions{i} = sprintf('Region %d from %s', i, obj.atlas_name);
        
    end
    
end

% Enforce row/column format (avoids compatibility problems later)

if iscolumn(obj.labels), obj.labels = obj.labels'; end

if ~iscolumn(obj.label_descriptions), obj.label_descriptions = obj.label_descriptions'; end
    
    

end