function atlas_obj = merge_atlases(atlas_obj, atlas_obj_to_add, varargin)
% Add all or some regions from an atlas object to another atlas object (with/without replacing existing labeled voxels)
%
%  atlas_obj = merge_atlases(atlas_obj, atlas_obj_to_add, [optional arguments])
%
% - First input atlas_obj defines space and voxel size
% - By default, values in atlas_obj_to_add are added to probability maps, if available,
%   and highest-probability region determines voxel label.  If probability maps are missing,
%   default is to replace labels in atlas_obj with the ones from the
%   to-be-added atlas.  The 'noreplace' option zeroes out the values in the to-add atlas if they are already
%   in the original atlas, privileging the atlas entered first.
%
%
% Optional arguments:
% Inputs to select_atlas_subset, e.g., vector of region numbers or cell array of region names
%  e.g., [1 3], {'VPL' 'IFG'}
%
% 'noreplace' : Do not replace voxels in original atlas; see above.
% 'always_replace' : always replace voxels in original atlas; see above.

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

doreplacevoxels = true;
doverbose = true;
always_replace = false;

% Find optional text strings and replace
for i = 1:length(varargin)
    
    if iscell(varargin{i})
        % will be passed in to select_atlas_subset
        
    elseif isnumeric(varargin{i})
        % will be passed in to select_atlas_subset
        
    elseif ischar(varargin{i})
        switch varargin{i}
            
            case 'always_replace', always_replace = true; varargin{i} = [];
            case 'noreplace', doreplacevoxels = false; varargin{i} = [];
            case 'noverbose', doverbose = false; varargin{i} = [];
                
            otherwise , warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% Subset
% -------------------------------------------------------------------------

if ~isempty(varargin) && any(cellfun(@isvector, varargin) | cellfun(@iscell, varargin))
    
    if doverbose, disp('Getting atlas subset with select_atlas_subset'); end
    
    args = varargin(cellfun(@isvector, varargin) | cellfun(@iscell, varargin));
    atlas_obj_to_add = select_atlas_subset(atlas_obj_to_add, args{:});
    
end

%% populate descriptive fields if mismatched
fname = {'labels', 'label_descriptions', 'labels_2', 'labels_3', 'labels_4', 'labels_5'};
n_regions = num_regions(atlas_obj);
n_regions_to_add = num_regions(atlas_obj_to_add);
for i = 1:length(fname)
    if any(~cellfun(@isempty, atlas_obj.(fname{i}))) || any(~cellfun(@isempty, atlas_obj_to_add.(fname{i}))) 
        if length(atlas_obj.(fname{i})) == n_regions || length(atlas_obj_to_add.(fname{i})) == n_regions_to_add
            if length(atlas_obj.(fname{i})) ~= n_regions
                atlas_obj.(fname{i}) = repmat({''}, 1, n_regions);
            end
            if length(atlas_obj_to_add.(fname{i})) ~= n_regions_to_add
                atlas_obj_to_add.(fname{i}) = repmat({''}, 1, n_regions_to_add);
            end
        end
    end
end
atlas_obj.label_descriptions = atlas_obj.label_descriptions(:);
atlas_obj_to_add.label_descriptions = atlas_obj_to_add.label_descriptions(:);
        

% -------------------------------------------------------------------------
% Resample space
% -------------------------------------------------------------------------
if doverbose, disp('Resampling space of atlas to add'); end

isdiff = compare_space(atlas_obj, atlas_obj_to_add);
if isdiff
    
    atlas_obj_to_add = resample_space(atlas_obj_to_add, atlas_obj);
    
end

% Make sure we have same voxel list
isdiff = compare_space(atlas_obj, atlas_obj_to_add);
if isdiff == 3
    atlas_obj = replace_empty(atlas_obj);
    atlas_obj_to_add = replace_empty(atlas_obj_to_add);
end

    
% -------------------------------------------------------------------------
% Add to atlas
% -------------------------------------------------------------------------

if doverbose
    disp('Adding regions'); 

    if ~doreplacevoxels
        disp('Not replacing voxels already in original atlas'); 
    end
end


n_regions = num_regions(atlas_obj);

has_pmaps = ~isempty(atlas_obj.probability_maps) && size(atlas_obj.probability_maps, 2) == n_regions;

toadd_n_regions = num_regions(atlas_obj_to_add);

toadd_has_pmaps = ~isempty(atlas_obj_to_add.probability_maps) && size(atlas_obj_to_add.probability_maps, 2) == toadd_n_regions;

if has_pmaps && toadd_has_pmaps
    % We have probability maps for both
    
    if always_replace    % eliminate voxels in original atlas that exist in new atlas to add
        
        wh = any(atlas_obj_to_add.probability_maps, 2); 
        atlas_obj.probability_maps(wh, :) = 0;
        
    end
    
    if ~doreplacevoxels  % eliminate voxels already in atlas
        
        wh = any(atlas_obj.probability_maps, 2);
        atlas_obj_to_add.probability_maps(wh, :) = 0;
    end
    
    % concatenating sparse matrices first is much faster than converting to full and the
    % concatenating full matrices so instead of doing a double conversion
    % and concatenating, try concatenating sparse matricies if either of
    % them are sparse and then converting to full
    %atlas_obj.probability_maps = [full(atlas_obj.probability_maps) full(atlas_obj_to_add.probability_maps)];
    if issparse(atlas_obj.probability_maps) || issparse(atlas_obj_to_add.probability_maps)
        if ~issparse(atlas_obj.probability_maps)
            atlas_obj.probability_maps = sparse(double(atlas_obj.probability_maps));
        end
        if ~issparse(atlas_obj_to_add.probability_maps)
            atlas_obj_to_add.probability_maps = sparse(double(atlas_obj_to_add.probability_maps));
        end
    end
    atlas_obj.probability_maps = full([atlas_obj.probability_maps atlas_obj_to_add.probability_maps]);
    
    atlas_obj = probability_maps_to_region_index(atlas_obj);
    
elseif has_pmaps && ~toadd_has_pmaps
    % Probability maps for first atlas only. Use these, preserving p maps
    
    if always_replace    % eliminate voxels in original atlas that exist in new atlas to add
        
        wh = any(atlas_obj_to_add.dat, 2);
        atlas_obj.probability_maps(wh, :) = 0;
        
    end
    
    % Reconstruct prob maps to add from index
    [p_recon, region_values] = condf2indic(atlas_obj_to_add.dat);
    if region_values(1) == 0, p_recon = p_recon(:, 2:end); end      % remove leading col for zeroes if there is one
    
    if ~doreplacevoxels  % eliminate voxels already in atlas
        
        wh = any(atlas_obj.probability_maps, 2);
        p_recon(wh, :) = 0;
        
    end
    
    atlas_obj.probability_maps = [atlas_obj.probability_maps double(p_recon)];
    
    atlas_obj = probability_maps_to_region_index(atlas_obj);
    
else  % if ~has_pmaps && ~toadd_has_pmaps
    % Probability maps for neither
    % Or Probability maps for to-add only; use index values in .dat either way
    
    wh = logical(atlas_obj_to_add.dat);  % voxels to replace
    
    if always_replace   % eliminate voxels in original atlas that exist in new atlas to add
        
        atlas_obj.dat(wh, :) = 0;
        
    end
    
    if ~doreplacevoxels  % eliminate voxels already in atlas
        wh(logical(atlas_obj.dat)) = false;
    end
    
    atlas_obj.dat(wh) = atlas_obj_to_add.dat(wh) + n_regions; % shift values to add to end of existing ones
    
end % cases

label_fields = {'labels' 'labels_2' 'labels_3' 'labels_4' 'labels_5'};
for f = label_fields
    fld = f{1};
    atlas_obj.(fld) = [atlas_obj.(fld) atlas_obj_to_add.(fld)];
end

atlas_obj.label_descriptions = [atlas_obj.label_descriptions(:); atlas_obj_to_add.label_descriptions(:)];

try
    atlas_obj.references = strvcat(atlas_obj.references, atlas_obj_to_add.references);
catch
    warning('Atlas references are not the same datatype. Attempting automatic conversion to char arrays. Please check resultant reference field.')
    atlas_obj.references = strvcat(char(atlas_obj.references), char(atlas_obj_to_add.references));
end


end % function

