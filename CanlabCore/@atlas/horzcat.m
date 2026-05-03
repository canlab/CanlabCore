function obj = horzcat(obj, varargin)
% horzcat Implements the horzcat ([a b]) operator on atlas objects.
%
% Combine two or more atlas-class objects by sequentially calling
% merge_atlases on each input. Equivalent to writing
% [atlas1 atlas2 atlas3].
%
% :Usage:
% ::
%
%     obj = horzcat(atlas1, atlas2, atlas3, ...)
%     obj = [atlas1 atlas2 atlas3 ...]    % equivalent
%
% :Inputs:
%
%   **obj:**
%        First atlas-class object; defines the space and voxel size.
%
%   **varargin:**
%        Additional atlas-class objects to merge into obj.
%
% :Outputs:
%
%   **obj:**
%        Atlas-class object containing the union of regions across
%        inputs (combined via merge_atlases).
%
% :See also:
%   - merge_atlases

for i = 1:length(varargin)
    
    obj = merge_atlases(obj, varargin{i});
    
end

end

