function atlas_obj = merge_atlases(atlas_obj, atlas_obj_to_add, varargin)
% Add all or some regions from an atlas object to another atlas object
%
%  atlas_obj = add_to_atlas(atlas_obj, atlas_obj_to_add, [optional arguments])
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

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

doreplacevoxels = true;
doverbose = true;

% Find optional text strings and replace
for i = 1:length(varargin)
    
    if iscell(varargin{i})
        % will be passed in to select_atlas_subset
        
    elseif isnumeric(varargin{i})
        % will be passed in to select_atlas_subset
        
    elseif ischar(varargin{i})
        switch varargin{i}
            
            case 'noreplace', doreplacevoxels = false; varargin{i} = [];
                
            otherwise , warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% Subset
% -------------------------------------------------------------------------

if ~isempty(varargin) && any(cellfun(@isvector, varargin) || cellfun(@iscell, varargin))
    
    if doverbose, disp('Getting atlas subset with select_atlas_subset'); end
    
    atlas_obj_to_add = select_atlas_subset(atlas_obj_to_add, varargin{:});
    
end

% -------------------------------------------------------------------------
% Resample space
% -------------------------------------------------------------------------
if doverbose, disp('Resampling space of atlas to add'); end

atlas_obj_to_add = resample_space(atlas_obj_to_add, atlas_obj);


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

toadd_has_pmaps = ~isempty(atlas_obj.probability_maps) && size(atlas_obj.probability_maps, 2) == toadd_n_regions;

% % % doreplacevoxels = false would be integrated here...

if has_pmaps && toadd_has_pmaps
    % We have probability maps for both
    
    if ~doreplacevoxels  % eliminate voxels already in atlas
        wh = any(atlas_obj.probability_maps, 2);
        atlas_obj_to_add.probability_maps(wh, :) = 0;
    end
    
    atlas_obj.probability_maps = [atlas_obj.probability_maps atlas_obj_to_add.probability_maps];
    
    atlas_obj = probability_maps_to_region_index(atlas_obj);
    
elseif has_pmaps && ~toadd_has_pmaps
    % Probability maps for first atlas only. Use these, preserving p maps
    
    % Reconstruct prob maps from index
    [p_recon, region_values] = condf2indic(atlas_obj_to_add.dat);
    if region_values(1) == 0, p_recon = p_recon(:, 2:end); end      % remove leading col for zeroes if there is one
    
    if ~doreplacevoxels  % eliminate voxels already in atlas
        wh = any(atlas_obj.probability_maps, 2);
        p_recon(wh, :) = 0;
    end
    
    atlas_obj.probability_maps = [atlas_obj.probability_maps p_recon];
    
    atlas_obj = probability_maps_to_region_index(atlas_obj);
    
else  % if ~has_pmaps && ~toadd_has_pmaps
    % Probability maps for neither
    % Or Probability maps for to-add only; use index values in .dat either way
    
    wh = logical(atlas_obj_to_add.dat);  % voxels to replace
    
    if ~doreplacevoxels  % eliminate voxels already in atlas
        wh(logical(atlas_obj.dat)) = false;
    end
    
    atlas_obj.dat(wh) = atlas_obj_to_add.dat(wh) + n_regions; % shift values to add to end of existing ones
    
end % cases

atlas_obj.labels = [atlas_obj.labels atlas_obj_to_add.labels];

atlas_obj.label_descriptions = [atlas_obj.label_descriptions; atlas_obj_to_add.label_descriptions];

atlas_obj.references = strvcat(atlas_obj.references, atlas_obj_to_add.references);


end % function

