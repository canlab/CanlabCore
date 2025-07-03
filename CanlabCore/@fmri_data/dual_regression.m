function [spatial_maps, timecourses, tmaps] = dual_regression(group_map_obj, data_obj, varargin)
% DUAL_REGRESSION Dual regression method for fmri_data objects
%
% :Usage:
% ::
%   [spatial_maps, timecourses, tmaps] = dual_regression(group_map_obj, data_obj, [optional inputs])
%
% :Inputs:
%   **group_map_obj**   : fmri_data object (e.g., ICA maps; the method is called on this object)
%   **data_obj** : fmri_data object 
%
% :Optional Inputs:
%   See dual_regression_fsl for options
%
% :Outputs:
%   **spatial_maps** : fmri_data object of spatial maps
%   **timecourses**  : matrix of component × timepoint
%   **tmaps**        : z-scored version of spatial maps (fmri_data object)
%
% Examples:
% -------------------------------------------------------------------------
% group_map_obj = load_image_set('npsplus');
% group_map_obj = get_wh_image(GM, 2:7);
% group_map_obj.image_names
% data_obj = fmri_data(which('swrsub-sid001567_task-pinel_acq-s1p2_run-03_bold.nii.gz'));
% 
% data_obj = resample_space(T, group_map_obj);
% [spatial_maps, timecourses, tmaps] = dual_regression_fsl(group_map_obj.dat, data_obj.dat);

data_obj = resample_space(data_obj, group_map_obj);

[spatial_maps, timecourses, tmaps] = dual_regression_fsl(group_map_obj, data_obj, varargin{:});


end