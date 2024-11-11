function [obj_subset, to_extract] = select_atlas_subset(obj, varargin)
% Select a subset of regions in an atlas by name or integer code, with or without collapsing regions together
%
% Note, some atlases are probablistic, and this may result in counterintuitive bheavior. If visualizing
% a probablistic atlas you only see p(region) > p(~region), a winner takes all parcellation scheme where
% a voxel is labeled as belonging to a region only if the probability is greater than the probability of
% it belonging to any other region. However, when you extract such a region you additionally get voxels
% with non-zero probability of belonging to said region even if other region labels are more probable.
% Consequently the extracted region boundaries will EXCEED those of the region you might have seen
% when visualizing the complete atlas, which can be confusing. Intelligent use of the probabilistic 
% labels is ideal, but if you prefer a more intuitive (albeit misleading) solution you should use the
% 'deterministic' option'
%
% Select a subset of regions in an atlas using:
% - a cell array of one or more strings to search for in labels
% - a vector of one or more integers indexing regions
%
% Output:
%       a new object, obj_subset
%       to_extract, a logical vector indicating which regions were selected
%
% obj_subset = select_atlas_subset(obj, varargin)
%
% Options:
% 'flatten' : Flatten integer vector in .dat to a single 1/0 mask. Good for
% creating masks from combinations of regions, or adding a set to another
% atlas as a single atlas region. .probability_maps is reestimated as
% posterior probability of union of constituent regions under conditional 
% independence assumption. If probability maps don't sum to 1 across voxels
% then conditional independence assumption is violated and we default to
% using the maximum probability instead.
%
% 'conditionally_ind' : If specified then 'flatten' assumes conditional
% independence even when probabilities don't sum to 1. Does nothing
% otherwise.
%
% 'labels_2' : If you enter any field name in the object, e.g., labels_2,
% the function will search for keyword matches here instead of in obj.labels.
%
% 'exact' : If you enter 'exact', function will not look for overlaps in
% names, and only look for exact string matches.
%
% 'regexp' : If you enter 'doregexp', function will treat your string as a
% regular expression.
%
% 'deterministic' : returns a labeled voxel iff p(region) > p(~region).
% Has no effect unless atlas has its probability_maps property populated,
% in which case the default behavior is to return all voxels with 
% p(region) > 0. Note: you probably want to apply a threshold operation too
% to remove low probability tissue boundary regions.
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
% Load a canonical atlas and select only the default mode, limbic, and brainstem regions
% based on the .labels_2 property
% -------------------------------------------------------
% atlas_obj = load_atlas('canlab2018_2mm');
% [obj_subset, to_extract] = select_atlas_subset(atlas_obj, {'Cerebellum', 'Def'}, 'labels_2');
% montage(obj_subset);
%
% see also: image_vector.select_voxels_by_value

% Programmers' notes:
% Created by tor wager
% Edited by tor, 12/2019, to use any field to select labels; and update
% auxiliary label fields too.
%
% Edited by Michael Sun, 09/04/2023, added an 'exact' flag, so that users
% can select to match string labels exactly and not capture regions with
% overlapping labels e.g., Cblm_Vermis_VI and Cblm_Vermis_VII.

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

strings_to_find = [];
integers_to_find = [];
doflatten = false;
condInd = false;
doexact=false;
doregexp=false;
deterministic=false;

% for entry of optional field to search for keywords
myfields = fieldnames(obj);
mylabelsfield = 'labels';

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
            case 'exact', doexact = true;

            case 'regexp', doregexp = true;

            case 'conditionally_ind', condInd = true;

            case 'deterministic', deterministic = true;
                
            otherwise
                
                if any(strcmp(varargin{i}, myfields))
                    mylabelsfield = varargin{i};
                    
                else
                    
                    warning(['Unknown input string option:' varargin{i} '. Assuming it might be an atlas label. Place atlas labels in a cell array']);
                    strings_to_find{1} = varargin{i};
                    
                end
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
    if doregexp == true
        wh = ~cellfun(@isempty, regexp(obj.(mylabelsfield), strings_to_find{i}));
        if strmatch(mylabelsfield, 'label_descriptions')
            wh=wh'; 
        end
    
        to_extract = to_extract | wh;

    elseif doexact == false
        wh = ~cellfun(@isempty, strfind(obj.(mylabelsfield), strings_to_find{i}));
        if strmatch(mylabelsfield, 'label_descriptions')
            wh=wh'; 
        end
    
        to_extract = to_extract | wh;
    % Addition by Michael Sun 09/04/2023: Match strings exactly
    else
        wh = strcmp(obj.(mylabelsfield), strings_to_find{i});
        if strmatch(mylabelsfield, 'label_descriptions')
            wh=wh'; 
        end
    
        to_extract = to_extract | wh;
    end
end

% -------------------------------------------------------------------------
% Check for problesum ms
% -------------------------------------------------------------------------

if ~isempty(obj.probability_maps) && size(obj.probability_maps, 2) ~= k  

    warning('The obj.probability_maps field is invalid because it has the wrong number of entries. Removing it from the created subset atlas.');
    obj.probability_maps = [];
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

% obj_subset.label_descriptions = obj_subset.label_descriptions(to_extract);

% if ~isempty(obj_subset.image_names) && size(obj_subset.image_names, 1) == k
%     obj_subset.image_names = obj_subset.image_names(to_extract, :);
% end

my_strings = {'image_names' 'fullpath'};

for i = 1:length(my_strings)
    
    if  ~isempty(obj_subset.(my_strings{i})) && size(obj_subset.(my_strings{i}), 1) == k
        obj_subset.(my_strings{i}) = obj_subset.(my_strings{i})(to_extract, :);
    end
    
end

my_strings = {'label_descriptions' 'image_names' 'fullpath' 'labels_2' 'labels_3' 'labels_4' 'labels_5'};

for i = 1:length(my_strings)
    if  ~isempty(obj_subset.(my_strings{i})) && length(obj_subset.(my_strings{i})) == k
        obj_subset.(my_strings{i}) = obj_subset.(my_strings{i})(to_extract);
    end
    
end

if ~isempty(obj.probability_maps) && size(obj.probability_maps, 2) == k  % valid p maps

    if deterministic
        % zero values where p(region) < p(~region)
        for this_region = find(to_extract)
            to_zero = any(obj_subset.probability_maps(:,this_region) < obj_subset.probability_maps(:, ~to_extract),2);
            obj_subset.probability_maps(to_zero,this_region) = 0;
        end
    end

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
    
    % this is wrong. Consider, you could have a voxel that has 50% 
    % probability of being voxel A and 50% probability of being voxel B.
    % the max value is %50, but the probability that it's one or the other
    % must be 100%. Let's check if it's safe to assume conditional 
    % independence instead.
    %
    % We don't know if the atlas labels are exhaustive, but if the sum of
    % atlas label probabilities exceeds 1 we can conclude that P(A,B) ~=
    % P(A) + P(B), because the total density must sum to 1 to be a valid
    % probability measure. We check if probabilities exceed 1 with some
    % tolerance to deal with floating point errors.
    hasdata = any(obj.probability_maps,2);
    tol = 1e-7;
    if ~condInd
        if all(sum(obj.probability_maps(hasdata,:), 2) - 1 < tol)
            condInd = true;
        end
    end
    if condInd
        obj_subset.probability_maps = sum(obj_subset.probability_maps, 2);
    else
        % we have evidence that probabilities are not conditionally
        % independent
        warning('Probabilities don''t sum to 1, so cannot assume condition independence. New map will use maximum probability instead of sum of probabilities.')
        obj_subset.probability_maps = max(obj_subset.probability_maps, [], 2);
    end

    if any(sum(obj_subset.probability_maps(hasdata,:), 2) - 1 > tol)
        warning('Max subset cumulative probability = %0.2f. Conditional independence assumption was probably violated.', max(sum(obj_subset.probability_maps(hasdata,:), 2)));
    end
    
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


