function obj = select_one_image(obj, wh)
% select_one_image Select a single image from a multi-image statistic_image object.
%
% Selects one image out of a set of images stored in a statistic_image
% object. Reduces .p, .ste, .sig, and .dat to a single column and resets
% .removed_images to 0.
%
% :Usage:
% ::
%
%     obj = select_one_image(obj, wh)
%
% :Inputs:
%
%   **obj:**
%        A statistic_image object.
%
%   **wh:**
%        An integer for which image in a series of images stored in obj
%        you want.
%
% :Outputs:
%
%   **obj:**
%        The input object reduced to the single image indexed by wh,
%        across the .p, .ste, .sig, and .dat fields.
%
% :See also:
%   - statistic_image
%   - get_wh_image

myfields = {'p' 'ste' 'sig' 'dat'};

for i = 1:length(myfields)
    
    if size(obj.(myfields{i}), 2) < wh
        
        error(['Field ' obj.(myfields{i}) ' does not have enough columns...object is incorrectly formed.']);
        
    end
    
    obj.(myfields{i}) = obj.(myfields{i})(:, wh);
    
end

%if obj.removed_images

obj.removed_images = 0;

%end


end
