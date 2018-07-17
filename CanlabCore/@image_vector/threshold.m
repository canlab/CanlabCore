function obj = threshold(obj, input_threshold, thresh_type, varargin)
% Threshold image_vector (or fmri_data or fmri_obj_image) object based on raw threshold values. 
% For statistical thresholding, convert to a statistic_image object and see the threshold method for that object.
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
%        Vector of 2 values defining data value bounds at which
%        to threshold, e.g., [0 Inf] or [-3 3]
%
%   **thresh_type:**
%        String: 'raw-between' or 'raw-outside'
%
% :Optional Inputs: Argument or argument followed by value:
%
%   **k:**
%        Followed by extent threshold cluster size, default = 1
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
%    % Retain positive values, cluster extent > 100 voxels
%    obj = threshold(obj, [0 Inf], 'raw-between', 'k', 100)
%
%    % Retain voxels with absolute value > 3
%    obj = threshold(obj, [-3 3], 'raw-outside')
%
% :See also:
% statistic_image.threshold, statistic_image.multi_threshold
%
% ..
%    Programmers' notes:
%    Tor: Updated documentation, July 2015
% ..

k = 1; % Inputs
dotrim = 0;
doverbose = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'k'
                k = varargin{i + 1};
                
                if isempty(obj.volInfo)
                    error('You must add a volInfo structure to the statistic image object to do extent-based thresholding');
                end
                
            case 'trim_mask'
                dotrim = 1;
                
            case 'noverbose'
                doverbose = 0;
                
            otherwise
                error('Illegal argument keyword.');
        end
    end
end

switch lower(thresh_type)
    
    case {'raw-between', 'raw-outside'}
        thresh = input_threshold;
        
    otherwise
        error('Unknown threshold type. Enter raw-between or raw-outside')
end

% -------------------------------------------------------
obj = replace_empty(obj);

switch lower(thresh_type)
    
    case 'raw-between'
        wh = obj.dat > thresh(1) & obj.dat < thresh(2);
        
        if doverbose
            fprintf('Keeping vals between %3.3f and %3.3f: %3.0f elements remain\n', thresh(1), thresh(2), sum(wh(:)));
        end
        
    case 'raw-outside'
        wh = obj.dat < thresh(1) | obj.dat > thresh(2);
        
        if doverbose
            fprintf('Keeping vals outside of %3.3f to %3.3f: %3.0f elements remain\n', thresh(1), thresh(2), sum(wh(:)));
        end
        
end

obj.dat(~wh) = 0;


% Apply size threshold
% --------------------------------------
if k > 1
    if doverbose, fprintf('Applying cluster extent threshold of %3.0f contiguous voxels\n', k); end
    for i = 1:size(obj.dat, 2)
        
        whkeep = logical(iimg_cluster_extent(double(obj.dat(:, i) ~= 0), obj.volInfo, k));
        obj.dat(~whkeep, i) = 0;
        
    end
end

obj = remove_empty(obj);

% Clean up and trim
% --------------------------------------

if dotrim
    obj = trim_mask(obj);
end


end % function

