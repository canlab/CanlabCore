function isdiff = compare_space(obj, obj2)
% Compare spaces of two image_vector objects
%
% :Usage:
% ::
%
%     function isdiff = compare_space(obj, obj2)
%
% Returns 0 if same, 1 if different spaces, 2 if no volInfo info for one or
% more objects. 3 if same space, but different in-mask voxels in .dat or
% volInfo.image_indx

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
