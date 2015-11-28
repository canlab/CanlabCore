function obj = trim_mask(obj)
% Exclude empty voxels from mask information in obj.volInfo structure, and re-make obj.volInfo
%
% ..
%    Tor Wager, 2013
% ..

obj = replace_empty(obj);
whomit = all(obj.dat == 0 | isnan(obj.dat), 2);

obj.dat(whomit, :) = [];
if ~isempty(obj.removed_voxels) && length(obj.removed_voxels) == obj.volInfo.n_inmask
    obj.removed_voxels(whomit) = [];
end
obj.volInfo.image_indx(obj.volInfo.wh_inmask(whomit)) = 0;
obj.volInfo.wh_inmask(whomit) = [];
obj.volInfo.xyzlist(whomit, :) = [];
obj.volInfo.cluster(whomit) = [];
obj.volInfo.n_inmask = length(obj.volInfo.wh_inmask);

if isa(obj, 'fmri_data')
    obj.mask.dat(whomit, :) = [];
    if ~isempty(obj.mask.removed_voxels) && length(obj.mask.removed_voxels) == obj.mask.volInfo.n_inmask
        obj.mask.removed_voxels(whomit) = [];
    end
    obj.mask.volInfo.image_indx(obj.mask.volInfo.wh_inmask(whomit)) = 0;
    obj.mask.volInfo.wh_inmask(whomit) = [];
    obj.mask.volInfo.xyzlist(whomit, :) = [];
    obj.mask.volInfo.cluster(whomit) = [];
    obj.mask.volInfo.n_inmask = length(obj.mask.volInfo.wh_inmask);
end

end % function
