function obj = threshold(obj, input_threshold, varargin)
% Threshold atlas object based on values in obj.probability_maps property
%
% - probability_maps property should contain values between 0 and 1.
% - This method will not work for empty probability_maps
%
% See image_vector.threshold and statistic_image.threshold for comparable
% methods for other object types. These do slightly different things,
% depending on object type.
%
% :Usage:
% ::
%
%    obj = threshold(obj, input_threshold, thresh_type, [optional arguments])
%
% This is a method for an image_vector object
%
% Thresholding is not reversible. For statistic_image objects it is.
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2015 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **obj:**
%        image_vector object
%
%   **input_threshold:**
%        Value between 0 and 1 defining probability threshold
%        OR
%        if thresh_type is 'range', Vector of 2 values defining
%        probability values to save
%        to threshold, e.g., [0 Inf] or [-3 3]
%
%   **thresh_type:**
%        Followed by one of the following strings:
%        'above' : keep prob. values above threshold; [default if no value entered]
%        'below' : keep prob. values below threshold
%        'range' : keep prob. values between input_threshold(1) and input_threshold(2)
%
% :Optional Inputs: Argument or argument followed by value:
%
%   **k:**
%        Followed by extent threshold cluster size, default = 1 (no cluster thresholding)
%
%   **remove_parcel_fragments:**
%        Probability maps can have multiple local maxima per parcel, resulting in fragmentation of 
%        parcels after thresholding. This flag eliminates fragments.
%
%   **spin_off_parcel_fragments:**
%        Probability maps can have multiple local maxima per parcel, resulting in fragmentation of 
%        parcels after thresholding. This flag creates new parcels from fragmented subparcels.
%
%   **trim_mask:**
%        Reduce the mask in obj.voInfo based on thresholding
%
%   **noverbose:**
%        Suppress verbose output
%
%
% :Outputs:
%
%   **obj:**
%        thresholded image_vector object
%
%
% :Examples:
% ::
%
% dat = threshold(dat, .2, 'k', 3);
% dat = threshold(dat, .2, 'k', 3, 'noverbose');
% dat = threshold(dat, .2, 'thresh_type', 'below', 'k', 3, 'noverbose');
%
%
% :See also:
% statistic_image.threshold, statistic_image.multi_threshold
%
% ..
%    Programmers' notes:
%    Tor: Updated documentation, July 2015
% ..

% Defaults and inputs
% --------------------------------------

doverbose = true;
dotrim = false;
k = 1;
thresh_type = 'above';
do_remove_frag = false;
do_spin_off_frag = false;

wh = strcmp(varargin, 'k');  
if any(wh), wh = find(wh); wh = wh(1); k = varargin{wh + 1}; varargin{wh + 1} = []; end
 
wh = strcmp(varargin, 'noverbose');  
if any(wh), doverbose = false; end

wh = strcmp(varargin, 'thresh_type');  
if any(wh), wh = find(wh); wh = wh(1); thresh_type = varargin{wh + 1}; varargin{wh + 1} = []; end

wh = strcmp(varargin,'remove_parcel_fragments');
if any(wh), do_remove_frag = true; end

wh = strcmp(varargin,'spin_off_parcel_fragments');
if any(wh), do_spin_off_frag = true;  end

% Error checking
% --------------------------------------

n = num_regions(obj); 

if isempty(obj.probability_maps) || size(obj.probability_maps, 2) ~= n
    error('Empty or invalid .probability_maps. Cannot threshold');
end

    
% Defaults and inputs
% --------------------------------------

switch thresh_type
    
    case 'above'
        
        if doverbose, fprintf('Keeping probability_map values above %3.2f\n', input_threshold); end
        
        obj.probability_maps(obj.probability_maps < input_threshold) = 0;
        
    case 'below'
        
        if doverbose, fprintf('Keeping probability_map values below %3.2f\n', input_threshold); end
        
        obj.probability_maps(obj.probability_maps > input_threshold) = 0;
        
    case 'range'
        
        if doverbose, fprintf('Keeping probability_map values between %3.2f and %3.2f\n', input_threshold); end
        
        wh = obj.probability_maps < input_threshold(1) || obj.probability_maps > input_threshold(2);
        obj.probability_maps(wh) = 0;
        
    otherwise
        
        error('Unrecognized thresh_type. Enter ''above'' ''below'' or ''range''')
        
end


% Re-compute parcel id vector
% --------------------------------------

obj = probability_maps_to_region_index(obj);

% Apply size threshold
% --------------------------------------
if k > 1
    
    if doverbose, fprintf('Applying cluster extent threshold of %3.0f contiguous voxels\n', k); end
    
    % iteratively remove clusters until removal operation has nothing left
    % to do
    obj = obj.replace_empty();
    newobj = apply_cluster_constraint(obj, k);
    while length(unique(newobj.dat)) ~= length(unique(obj.dat))
        obj = newobj;
        newobj = apply_cluster_constraint(obj, k);
    end
    obj = newobj;
end


if do_spin_off_frag || do_remove_frag
    obj = spin_off_frag(obj);
end

if do_remove_frag
    obj = remove_frag(obj);
end

% Clean up and trim
% --------------------------------------

if dotrim
    obj = trim_mask(obj);
end

end % function


% % use image vector threshold
% % this thresholds the dat field, so enter probability values into .dat
% obj_tmp = obj;
% obj_tmp.dat = obj_tmp.probability_maps;
% 
% obj_tmp = image_vector.threshold(obj_tmp, input_threshold, thresh_type, varargin{:});
% 
% for i = 1:k
% dat_wh = get_wh_image(dat, 1); 
% dat_wh = threshold(dat_wh, [.2 Inf], 'raw-between', 'k', 3);

function obj = apply_cluster_constraint(obj, k)
    new_obj = obj;
    new_obj.probability_maps = [];

    uniq_roi = unique(new_obj.dat, 'stable');
    uniq_roi(uniq_roi == 0) = [];

    for i = 1:length(uniq_roi)
        this_obj = new_obj.select_atlas_subset(i);
        
        anyvalue = any(double(this_obj.dat), 2);
            
        whkeep = logical(iimg_cluster_extent(anyvalue, obj.volInfo, k));
            
        obj.probability_maps(~whkeep, i) = 0;
    end

    obj = obj.probability_maps_to_region_index;
end


% Identifies any fragmented parcels, splits them in two and postfixes
% _fragN to the label name, retaining corresponding probabilities for
% each regoin. N is the fragment index, frag1 is the largest and subsequent
% parcels are in order of size.
function obj = spin_off_frag(obj)
    new_obj = obj;
    new_obj.probability_maps = [];

    uniq_roi = unique(new_obj.dat, 'stable');
    uniq_roi(uniq_roi == 0) = [];

    pmaps = [];
    labels = {};
    extra_labels = struct('labels_2',{{''}}, 'labels_3',{{''}}, 'labels_4', {{''}}, 'labels_5', {{''}}, 'label_descriptions', {{''}});
    for i = 1:length(uniq_roi)
        this_obj = new_obj.select_atlas_subset(i);
        
        r = atlas2region(this_obj);
        r = r.reparse_continguous();

        numVox = arrayfun(@(x1)(x1.numVox), r.reparse_continguous);
        [numVox, I] = sort(numVox,'descend');
        for j = 1:length(numVox)
            new_region = r(I(j)).region2atlas(obj);
            pmaps = [pmaps, zeros(size(obj.probability_maps,1),1)];
            pmaps(new_region.dat == 1,end) = obj.probability_maps(new_region.dat == 1,i);
            if length(numVox) > 1
                labels{end+1} = [obj.labels{i},'_frag',num2str(j)];
                fprintf('Spinning off %s\n', labels{end});
            else
                labels{end+1} = obj.labels{i};
            end
            for fname = fieldnames(extra_labels)'
                if length(obj.(fname{1})) == length(obj.labels)
                    if strcmp(extra_labels.(fname{1}){1},'')
                        extra_labels.(fname{1}) = obj.(fname{1})(i);
                    else
                        extra_labels.(fname{1}){end+1} = obj.(fname{1}){i};
                    end
                end
            end
        end
    end

    obj.probability_maps = pmaps;
    obj.labels = labels;
    obj = obj.probability_maps_to_region_index;
    for fname = fieldnames(extra_labels)'
        obj.(fname{1}) = extra_labels.(fname{1});
    end

    obj.label_descriptions = obj.label_descriptions(:);
end

% removes all parcels after the largest parcel fragment
function obj = remove_frag(obj)
    isfrag = contains(obj.labels, '_frag');
    isPrincipalFrag = contains(obj.labels,'_frag1');
    keep = isPrincipalFrag | ~isfrag;
    
    if sum(~keep) > 0
        fprintf('Removing %s\n', obj.labels{~keep});
    end
    obj = obj.select_atlas_subset(find(keep));
    obj.labels = cellfun(@(x1)strrep(x1,'_frag1',''), obj.labels, 'UniformOutput', false);
end