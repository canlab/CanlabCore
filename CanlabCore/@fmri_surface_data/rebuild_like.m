function newobj = rebuild_like(obj, newdat)
% rebuild_like Build a new fmri_surface_data carrying obj's geometry, with new data.
%
% :Usage:
% ::
%     newobj = rebuild_like(obj, newdat)
%
% Helper used by data-transforming methods (mean, ica, predict, threshold, ...)
% to wrap a new [grayordinates x maps] data matrix back into an
% fmri_surface_data that carries the original brain_model / geom / intent /
% surface_space -- the surface analogue of fmri_data's volInfo-based re-wrap.
% This avoids the inherited image_vector/fmri_data rebuild that would emit a
% volume object.
%
% The new data must have the same number of grayordinate rows as obj (the
% geometry is unchanged); the number of maps (columns) may differ. Per-map
% annotations (image_names/X/Y/...) are kept only if the map count is unchanged.
%
% :Inputs:
%   **obj:**    template fmri_surface_data object (geometry source).
%   **newdat:** [nGrayordinates x K] numeric data.
%
% :Outputs:
%   **newobj:** fmri_surface_data with .dat = single(newdat), same geometry.
%
% :See also: fmri_surface_data, mean, reconstruct_image

if size(newdat, 1) ~= size(obj.dat, 1)
    error('fmri_surface_data:rebuild_like:rowmismatch', ...
        ['newdat has %d rows but the object has %d grayordinates. rebuild_like ' ...
         'preserves geometry; row count must match.'], size(newdat,1), size(obj.dat,1));
end

newobj = obj;
newobj.dat = single(newdat);

K = size(newdat, 2);
newobj.removed_voxels = false(size(newdat, 1), 1);
newobj.removed_images = false(K, 1);

% Drop per-map annotations if the number of maps changed (they no longer align)
if K ~= size(obj.dat, 2)
    newobj.image_names = {};
    newobj.X = [];
    newobj.Y = [];
    newobj.covariates = [];
    newobj.metadata_table = [];
end
end
