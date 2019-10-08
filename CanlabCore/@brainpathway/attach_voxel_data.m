% This is a stub with an example of how to attach voxel data, which will
% trigger the calculation of derivatives (region averages)
% This function will likely not be necessary, as the functions can be accomplished
% with listeners - this is an example really.

% img = fmri_data('/Users/torwager/Downloads/smoothed_residual.nii.gz')
% b = brainpathway();

% Resample image space of data to the atlas space
% (We don't want to do it the other way around, because it will make it
% difficult to have a consistent atlas definition and aggregate across
% different datasets)

img_resampled = resample_space(img, b.region_atlas);

% Assign the resampled data to the brainpathway object
% This triggers several steps:
% Updating node and region averages, and node and region connectivity
b.voxel_dat = img_resampled.dat;

% Plot the inter-region connectivity
% plot_connectivity(b);

% Load a different atlas instead, with a different parcellation:
% This triggers updates of the regions and region connectivity
b.region_atlas = load_atlas('yeo17networks');

b = brainpathway(load_atlas('yeo17networks'));
b.voxel_dat = img_resampled.dat;
