function obj = attach_voxel_dat(obj, img_obj) %, src, evt)
% attach_voxel_dat Resample and attach voxel-level data to the brainpathway object.
%
% :Usage:
% ::
%     obj = obj.attach_voxel_dat(img_obj);  % needs work: , src, evt)
%
% :Inputs:
%
%   **obj:**
%        A brainpathway handle object.
%
%   **img_obj:**
%        An fmri_data object containing image data to be resampled.
%
%   **src:**
%        Source object triggering the event (unused in this implementation).
%
%   **evt:**
%        Event data (unused in this implementation).
%
% :Outputs:
%
%   **obj:**
%        The updated brainpathway object with the voxel_dat property set to the
%        resampled image data.
%
% :Description:
%
%   This method uses the function resample_space to resample the space of the input
%   fmri_data object (img_obj) to the atlas space defined in the brainpathway object's
%   region_atlas property. The resampled data (img_resampled.dat) is then assigned
%   to obj.voxel_dat, triggering further updates such as recalculation of node and
%   region averages and connectivity.
%
% :Examples:
% ::
% % load an atlas that we'd like to use
% atlas_obj = load_atlas('canlab2024');
%
% % Create a brainpathway object with this atlas 
% % (if we omit atlas_obj, it will use a default atlas instead).
% brainpathway_obj = brainpathway(atlas_obj);
%
% % Resample image data in the fmri_data obj_denoised and add it to the
% % .voxel_data field in the brainpathway object, triggering updates of the
% % extracted region averages and connectivity matrix 
% brainpathway_obj.attach_voxel_data(obj_denoised);
%
% % Plot the region x region connectivity matrix 
% plot_connectivity(brainpathway_obj, 'regions');
% plot_connectivity(brainpathway_obj, 'regions', 'partitions', brainpathway_obj.region_atlas.labels_3);
%
% See also: resample_space

% Resample the image space of the provided fmri_data object to the atlas space.
% (We don't want to do it the other way around, because it will make it
% difficult to have a consistent atlas definition and aggregate across
% different datasets)

% img is an fmri_data object with data to attach. the data are in the .dat
% field.  These should be assigned to obj.voxel_dat

img_resampled = resample_space(img_obj, obj.region_atlas);

% Update the voxel_dat property with the resampled data.
% This triggers several steps:
% Updating node and region averages, and node and region connectivity
obj.voxel_dat = img_resampled.dat;

end


% examples
% Load a different atlas instead, with a different parcellation:
% This triggers updates of the regions and region connectivity
% b.region_atlas = load_atlas('yeo17networks');
%
% b = brainpathway(load_atlas('yeo17networks'));
% b.voxel_dat = img_resampled.dat;
