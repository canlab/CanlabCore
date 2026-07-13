function proxy = as_fmri_data_proxy(obj)
% as_fmri_data_proxy Wrap a surface object as an fmri_data with a dummy 1-D volInfo.
%
% Lets geometry-agnostic fmri_data analysis methods (predict, ica, ttest, ...) be
% reused on grayordinate data: each grayordinate is treated as a "voxel" in a
% 1-D volume. The returned proxy carries the same .dat / .X / .Y / covariates so
% the analysis is identical; only spatial reconstruction (which the analysis does
% not need) is meaningless on the proxy. Map geometry-bearing results back to a
% surface object with rebuild_like.
%
% :Inputs:  **obj:** an fmri_surface_data object.
% :Outputs: **proxy:** an fmri_data object [nGrayordinates x nMaps].
%
% :See also: predict, ica, ttest, rebuild_like, fmri_surface_data

nGray = size(obj.dat, 1);

vi = struct('mat', eye(4), 'dim', [nGray 1 1], 'dt', [16 0], ...
    'xyzlist', [(1:nGray)' ones(nGray,1) ones(nGray,1)], 'nvox', nGray, ...
    'image_indx', true(nGray,1), 'wh_inmask', (1:nGray)', 'n_inmask', nGray, 'fname', '');

iv = image_vector;
iv.volInfo = vi;
iv.dat = double(obj.dat);
iv.removed_voxels = false(nGray, 1);
iv.removed_images = false(size(obj.dat, 2), 1);
iv.image_names = obj.image_names;

proxy = fmri_data(iv);

% Carry the per-map analysis annotations
if ~isempty(obj.Y),          proxy.Y = obj.Y; end
if ~isempty(obj.X),          try, proxy.X = obj.X; catch, end; end %#ok<NOCOM>
if ~isempty(obj.covariates), proxy.covariates = obj.covariates; end
if ~isempty(obj.Y_names),    proxy.Y_names = obj.Y_names; end
end
