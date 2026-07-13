function tf = isempty(obj)
% isempty True only when an fmri_surface_data has no grayordinate data.
%
% Overrides image_vector/isempty, which also treats an empty volInfo as "empty".
% A cortex-only surface object legitimately has an EMPTY volInfo (only the
% subcortical voxel sub-block populates volInfo; see design decision D3), so the
% inherited test wrongly reports a fully-populated cortical map as empty. Here an
% object is empty only when its .dat has no rows.
%
% :See also: fmri_surface_data, image_vector.isempty

tf = builtin('isempty', obj.dat) || size(obj.dat, 1) == 0;
end
