function obj = enforce_variable_types(obj)
% Re-casts variables in objects into standard data types, which can save space
% in memory.  Also removes empty voxels as a space-saving device.
%
% obj = enforce_variable_types(obj)
%
% Needs testing to ensure that doing this does not break any other
% display/computation functions.

obj = remove_empty(obj);
obj.dat = single(obj.dat);

obj.volInfo.wh_inmask = uint32(obj.volInfo.wh_inmask);  % 4 billion positive valued integers
obj.volInfo.xyzlist = uint16(obj.volInfo.xyzlist);
obj.volInfo.cluster = uint32(obj.volInfo.cluster);

if isa(obj, 'fmri_data')
    % future: remove .mask altogether...?
    obj.mask.dat = single(obj.mask.dat);
    
    obj.mask = remove_empty(obj.mask);
    
    obj.mask.volInfo.wh_inmask = uint32(obj.mask.volInfo.wh_inmask);  % 4 billion positive valued integers
    obj.mask.volInfo.xyzlist = uint16(obj.mask.volInfo.xyzlist);
    obj.mask.volInfo.cluster = uint32(obj.mask.volInfo.cluster);
    
end

if isa(obj, 'statistic_object')
    
    obj.sig = logical(obj.sig);
    
    obj.p = single(obj.p);
    
end

end
