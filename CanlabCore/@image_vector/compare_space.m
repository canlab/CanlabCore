function isdiff = compare_space(obj, obj2)
% compare_space Compare image spaces of two image_vector objects.
%
% Compares the affine matrix, dimensions, and voxel count of two
% image_vector / fmri_data objects, plus their in-mask voxel
% configuration.
%
% :Usage:
% ::
%
%     isdiff = compare_space(obj, obj2)
%
% :Inputs:
%
%   **obj:**
%        First image_vector / fmri_data object.
%
%   **obj2:**
%        Second image_vector / fmri_data object.
%
% :Outputs:
%
%   **isdiff:**
%        Integer code:
%
%        - 0 if same
%        - 1 if different spaces
%        - 2 if no volInfo info for one or more objects
%        - 3 if same space, but different in-mask voxels in .dat or
%          volInfo.image_indx
%
% :Examples:
% ::
%
%     d = compare_space(dat, mask);
%
% :See also:
%   - resample_space
%   - apply_mask

if isempty(obj.volInfo) || isempty(obj2.volInfo)
    isdiff = 2;
    return;
end

n1 = [obj.volInfo.mat(:); obj.volInfo.dim(:); obj.volInfo.nvox];

n2 = [obj2.volInfo.mat(:); obj2.volInfo.dim(:); obj2.volInfo.nvox];

isdiff = any(n1 - n2);

if ~isdiff
    if size(obj.dat, 1) ~= size(obj2.dat, 1) || any(obj.volInfo.image_indx - obj2.volInfo.image_indx)
       isdiff = 3; 
    end
    
end

end
