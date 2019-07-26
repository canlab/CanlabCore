function obj = horzcat(obj, varargin)
% Use merge_atlas to combine atlases, or newatlas = (atlas1, atlas2, atlas3);

for i = 1:length(varargin)
    
    obj = merge_atlases(obj, varargin{i});
    
end

end

