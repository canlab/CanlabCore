function obj = select_one_image(obj, wh)
%
% Selects one image out of a set of images stored in a statistic_image object
%
% :Inputs:
%
%   **obj:**
%        A statistic_image object
%
%   **wh:**
%        An integer for which image in a series of images stored in obj you want

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
