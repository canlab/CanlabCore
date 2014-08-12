function obj = threshold(obj, input_threshold, thresh_type, varargin)
%
% Threshold image_vector (or fmri_data or fmri_obj_image) object based on
% raw threshold values. For statistical thresholding, convert to a
% statistic_image object and see the threshold method for that object.
%
% Examples:
% -> Retain positive values, cluster extent > 100 voxels
% obj = threshold(obj, [0 Inf], 'raw-between', 'k', 100)
%
% -> Retain voxels with absolute value > 3
% obj = threshold(obj, [-3 3], 'raw-outside')

% Inputs
% -------------------------------------------------------
k = 1;
dotrim = 0;

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
        
        fprintf('Keeping vals between %3.3f and %3.3f: %3.0f elements remain\n', thresh(1), thresh(2), sum(wh(:)));
        
    case 'raw-outside'
        wh = obj.dat < thresh(1) | obj.dat > thresh(2);
        
        fprintf('Keeping vals outside of %3.3f to %3.3f: %3.0f elements remain\n', thresh(1), thresh(2), sum(wh(:)));
        
end

obj.dat(~wh) = 0;


% Apply size threshold
% --------------------------------------
if k > 1
    fprintf('Applying cluster extent threshold\n');
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

