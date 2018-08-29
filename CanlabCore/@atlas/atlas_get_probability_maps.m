function prob_map_obj = atlas_get_probability_maps(obj, varargin)
% Given an atlas object, return probability maps or indicator maps based on parcels
%
% :Usage:
% ::
% prob_map_obj = atlas_get_probability_maps(obj)
%
% :Inputs:
%
%   **obj:**
%        an atlas-class object
%
% :Outputs:
%
%   **prob_map_obj:**
%        an fmri_data object containing the probability maps or indicators
%        as data
% 
% Copyright 2018 Tor Wager

% Get maps if they exist
% -------------------------------------------------------------------------

if isempty(obj.probability_maps)
    
    obj.probability_maps = single(condf2indic(obj.dat, 'integers'));
    
    obj.probability_maps(isnan(obj.probability_maps)) = 0;
    
end

% Get maps and construct an object
% -------------------------------------------------------------------------

prob_map_obj = image_vector('dat', single(full(obj.probability_maps)), ...
    'volInfo', obj.volInfo, ...
    'removed_voxels', obj.removed_voxels, 'removed_images', obj.removed_images, ...
    'image_names', obj.image_names, 'noverbose');

prob_map_obj = fmri_data(prob_map_obj);


end % function

