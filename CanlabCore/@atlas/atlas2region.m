function r = atlas2region(atlas_obj, varargin)
% Convert an atlas object to a region object
% r = atlas2region(atlas_obj)
%
% Note: There are multiple potential ways this can work.
% The default is to use uses the atlas_obj.dat field, 
% which has one label per voxel. For some atlases, a region may have no
% voxels with the max probability label. 
% 
% 'use_probabilities'  uses atlas_obj.probability maps if available, 
% which guarantees at least one region per atlas region. If atlas.probability_maps is
% empty, however, the function defaults to using index labels.
% Note: this can be slow for large atlases.
%
% Another default is to parse regions with each label into contiguous
% regions. This can result in multiple contiguous regions in the region object (r)
% for a single atlas region. (e.g., those in left and right hemispheres, if
% they are given the same label in the atlas).
%
% 'nocontiguous' (or 'unique_mask_values') is an option that suppresses the 
% parsing, creating one region per label/atlas region.
% If atlas.probability_maps is empty or 'use_labels' is specified,
% the function defaults to 'nocontiguous'.
%
% Examples:
% r = atlas2region(atlas_obj);                       % Use indices in dat.dat, one region per atlas region.
% r = atlas2region(atlas_obj, 'use_probabilities');  % Use probability_maps if available  
% r = atlas2region(atlas_obj, 'nocontiguous');       % One region per atlas region
%
% Tor Wager, Jan 2018

% Defaults and inputs
% ----------------------------------------------------------------------------

use_probabilities_if_available = false;
region_parse_method = 'contiguous_regions';

for i = 1:length(varargin)
    
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'use_probabilities', use_probabilities_if_available = false;
            case {'nocontiguous', 'unique_mask_values'}
                region_parse_method = 'unique_mask_values';
                
            otherwise , warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Checks
% ----------------------------------------------------------------------------

n_regions = num_regions(atlas_obj);
has_pmaps = ~isempty(atlas_obj.probability_maps) && size(atlas_obj.probability_maps, 2) == n_regions;

% Parse into regions
% ----------------------------------------------------------------------------

if use_probabilities_if_available && has_pmaps
    % Use probability_maps
    % --------------------------------------------------------------------
    r = {};
    wh_label = {};
    
    for i = 1:n_regions
        
        atlasregion = select_atlas_subset(atlas_obj, i);
        
        % recast variable types: avoid conflicts, e.g., in spm orthviews
        atlasregion.dat = double(atlasregion.dat);
        
        % parse into contiguous or index based, depending on inputs
        r{i} = region(atlasregion, region_parse_method, 'noverbose');
        
        % assign label indices, depending on how many regions we have for
        % this atlas region.
        wh_label{i} = i * ones(1, length(r{i}));
        
    end
    
    % flatten into single region object
    r = cat(2, r{:}); 
    wh_label = cat(2, wh_label{:});
    
else
    % Use region labels in atlas_obj.dat
    % --------------------------------------------------------------------
    
    atlas_obj.dat = double(atlas_obj.dat); % int32 doesn't work - conflict when adding vars to orthviews of diff types
    
    r = region(atlas_obj, 'unique_mask_values', 'noverbose');
    
    wh_label = 1:n_regions;    % for label assignment
    
    if n_regions > length(atlas_obj.labels)
        error('More regions than labels! Labels or integer codes may be wrong.');
    end
    
end % Options, pmap or .dat based


% Add labels to region names
% ----------------------------------------------------------------------------

r(1).custom_info1_descrip = atlas_obj.space_description;

if isempty(atlas_obj.labels)
    return
    
else
    
    for i = 1:length(r)

        r(i).shorttitle = atlas_obj.labels{wh_label(i)};
        
        if ~isempty(atlas_obj.label_descriptions) && length(atlas_obj.label_descriptions) >= i
            
            r(i).title = sprintf('%s from %s', atlas_obj.label_descriptions{wh_label(i)}, atlas_obj.atlas_name);
            
        end
        
    end
    
end

end

