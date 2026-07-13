function obj = replace_empty(obj, varargin)
% replace_empty (fmri_surface_data): no-op. .dat is always full-size.
%
% :Usage:
% ::
%     obj = replace_empty(obj)
%
% For image_vector/fmri_data, replace_empty re-expands a space-reduced .dat back
% to its original size (re-inserting removed voxels/images as zeros) so the data
% can be reconstructed into a volume. fmri_surface_data never reduces .dat -- it
% always holds the full [nGrayordinates x nMaps] grayordinate set in 1:1 row
% correspondence with .brain_model -- so there is nothing to re-insert.
%
% This override returns the object UNCHANGED, which removes the entire class of
% "forgot to replace_empty before reconstructing" bugs. Methods such as
% reconstruct_image / write / cat therefore need no replace_empty pre-step.
%
% See docs/fmri_surface_data_design_plan.md (decision D5b).
%
% :See also: remove_empty, fmri_surface_data, reconstruct_image

% Intentionally a no-op.
end
