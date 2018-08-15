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

wh = strcmp(varargin, 'k');  
if any(wh), wh = find(wh); wh = wh(1); k = varargin{wh + 1}; varargin{wh + 1} = []; end
 
wh = strcmp(varargin, 'noverbose');  
if any(wh), doverbose = false; end

wh = strcmp(varargin, 'thresh_type');  
if any(wh), wh = find(wh); wh = wh(1); thresh_type = varargin{wh + 1}; varargin{wh + 1} = []; end


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


% Apply size threshold
% --------------------------------------
if k > 1
    
    if doverbose, fprintf('Applying cluster extent threshold of %3.0f contiguous voxels\n', k); end
    
    for i = 1:size(obj.dat, 2)
        
        anyvalue = any(double(obj.probability_maps), 2);
        
        whkeep = logical(iimg_cluster_extent(anyvalue, obj.volInfo, k));
        
        obj.probability_maps(~whkeep, :) = 0;
        
    end
end


% Re-compute parcel id vector
% --------------------------------------

obj = probability_maps_to_region_index(obj);


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
