function fmri_obj = to_fmri_data(obj)
% to_fmri_data Export the subcortical/volumetric grayordinates to an fmri_data object.
%
% :Usage:
% ::
%     fmri_obj = to_fmri_data(surf_obj)
%
% Returns the volumetric (subcortical / cerebellar) part of a grayordinate
% fmri_surface_data object as a standard CANlab fmri_data object in the CIFTI
% volume's MNI space. Surface (cortical) grayordinates have no voxel coordinates
% and are dropped here -- use surf2vol (M4) to project cortex into a volume.
%
% The resulting fmri_data can be written to a NIfTI (.nii) with its write
% method, montaged, resampled, etc. like any volumetric object.
%
% :Inputs:
%   **obj:** an fmri_surface_data object with a volumetric sub-block
%            (brain_model.vol + one or more 'vox' models).
%
% :Outputs:
%   **fmri_obj:** an fmri_data object [n_subcortical_voxels x nMaps].
%
% :Examples:
% ::
%     s = fmri_surface_data(which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii'));
%     vol = to_fmri_data(s);
%     % vol.write('sample_filename', '/tmp/subctx.nii');   % to disk
%
% :See also: fmri_surface_data, reconstruct_image, build_volinfo_subblock, surf2vol

[vi, rows] = build_volinfo_subblock(obj.brain_model);
if isempty(vi) || isempty(rows)
    error('fmri_surface_data:to_fmri_data:novolume', ...
        ['This object has no volumetric (subcortical) grayordinates to export. ' ...
         'Use surf2vol to project cortical surface data into a volume.']);
end

iv = image_vector;
iv.volInfo = vi;
iv.dat = single(obj.dat(rows, :));
iv.removed_voxels = false(size(iv.dat, 1), 1);
iv.removed_images = false(size(iv.dat, 2), 1);
iv.image_names = obj.image_names;
iv.history = obj.history;
iv.history{end+1} = sprintf('to_fmri_data: extracted %d subcortical voxels x %d maps from %s', ...
    size(iv.dat,1), size(iv.dat,2), class(obj));

fmri_obj = fmri_data(iv);
end
