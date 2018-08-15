function [obj_subset, to_extract] = select_atlas_subset(obj, varargin)
% Select a subset of regions in an atlas by name or integer code, with or without collapsing regions together
%
% Select a subset of regions in an atlas using:
% - a cell array of one or more strings to search for in labels
% - a vector of one or more integers indexing regions
% Return output in a new object, obj_subset
%
% obj_subset = select_atlas_subset(obj, varargin)
%
% Options:
% 'flatten' : Flatten integer vector in .dat to a single 1/0 mask. Good for
% creating masks from combinations of regions, or adding a set to another
% atlas as a single atlas region. .probability_maps reduced to single max
% map
%
% Examples:
% 
% atlasfile = which('Morel_thalamus_atlas_object.mat');
% load(atlasfile)
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, [1 3]);
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, [1 3], {'VPL'});
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'VPL'});             % Sensory VPL thalamus, when using morel atlas
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'VPL' 'VPM' 'VPI'}); % Sensory thalamus, both lat and medial
%
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'Hb'});              % Habenula
%
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'CL' 'CeM' 'CM' 'Pf'});  % Intralaminar
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'Pv' 'SPf'});           % Midline group
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'CL' 'CeM' 'CM' 'Pf' 'Pv' 'SPf'}); % Intralaminar and Midline group
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'MD'});              % Mediodorsal nuc.
%
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'Pu'});              % Pulvinar
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'LGN'});             % LGN
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'MGN'});             % MGN
%
% Create a mask from a set of regions:
% obj_subset = select_atlas_subset(atlas_obj, {'CL' 'CeM' 'CM' 'Pf'}, 'flatten');  % Intralaminar
% r = atlas2region(obj_subset);
% orthviews(r)
%
% see also: image_vector.select_voxels_by_value

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

strings_to_find = [];
integers_to_find = [];
doflatten = false;

% optional inputs with default values
for i = 1:length(varargin)
    
    if iscell(varargin{i})
        strings_to_find = varargin{i};
        
    elseif isnumeric(varargin{i})
        integers_to_find = varargin{i};
        
    elseif ischar(varargin{i})
        switch varargin{i}

             case 'flatten', doflatten = true;

                 %             case 'xxx', xxx = varargin{i+1}; varargin{i+1} = [];
                
            otherwise , warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% INIT
% -------------------------------------------------------------------------

k = num_regions(obj);
to_extract = false(1, k);

% -------------------------------------------------------------------------
% FIND BY STRING
% -------------------------------------------------------------------------

for i = 1:length(strings_to_find)
    
    % Find which names match
    wh = ~cellfun(@isempty, strfind(obj.labels, strings_to_find{i}));
    
    to_extract = to_extract | wh;
    
end

% -------------------------------------------------------------------------
% FIND BY NUMBERS
% -------------------------------------------------------------------------

to_extract(integers_to_find) = true;

if ~any(to_extract)
    error('No regions identified to extract.');
end

% -------------------------------------------------------------------------
% Extract
% -------------------------------------------------------------------------

obj_subset = obj;

% labels and names

obj_subset.labels = obj_subset.labels(to_extract);
obj_subset.label_descriptions = obj_subset.label_descriptions(to_extract);

if ~isempty(obj_subset.image_names) && size(obj_subset.image_names, 1) == k
    obj_subset.image_names = obj_subset.image_names(to_extract, :);
end

    
if ~isempty(obj.probability_maps) && size(obj.probability_maps, 2) == k  % valid p maps

    obj_subset.probability_maps(:, ~to_extract) = [];
    
    obj_subset = probability_maps_to_region_index(obj_subset);
    
else % must use .dat vector with integers
    
    % rebuild integer index
    
    dat = zeros(size(obj_subset.dat));
    wh_cols = find(to_extract);
    
    for i = 1:length(wh_cols)
        
        wh_vals = obj_subset.dat == wh_cols(i);
        
        dat(wh_vals) = i;
        
    end
    
    obj_subset.dat = dat;
    
end

if doflatten
    % Flatten into mask
    % -------------------------------------------------------------------
    obj_subset.dat = single(obj_subset.dat > 0);
    
    obj_subset.probability_maps = max(obj_subset.probability_maps, [], 2);
    
    % Add labels for combined mask
    
    mylabels = [];
    if iscell(strings_to_find)
        mylabels = cat(2, strings_to_find{:});
    end
    
    if ~isempty(integers_to_find)
        morelabels = cat(2, obj.labels{integers_to_find});
        mylabels = [mylabels morelabels];
    end
    
    obj_subset.labels = {mylabels};
    
    obj_subset.label_descriptions = {sprintf('Combined subset: %s', mylabels)};
    
        % Make sure they are length of prob maps, to conform to standard atlas
    % object
%     np = size(obj_subset.probability_maps, 2);
%     if np
%         obj_subset.labels = repmat(obj_subset.labels, 1, np);
%         obj_subset.label_descriptions = repmat(obj_subset.label_descriptions, np, 1);
%     end
    
end % flatten

end % function


