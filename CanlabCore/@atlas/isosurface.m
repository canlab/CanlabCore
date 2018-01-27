function surface_handles = isosurface(atlas_obj, varargin)
% Create a series of surfaces in different colors, one for each region
% - Options for single color
%
% surface_handles = isosurface(atlas_obj, [optional arguments])
%
% optional arguments: 
% - Any optional inputs to imageCluster, e.g., 'alpha'
% - 'colors', followed by single color in { } or cell array of multiple colors
% - 'nomatchleftright', do not match colors across hemispheres (left/right)
%               Note: The default matches, and may override your colors.
%
% Examples:
% atlasfile = which('Morel_thalamus_atlas_object.mat');
% load(atlasfile)
%
% surface_handles = isosurface(atlas_obj);
% surface_handles = isosurface(atlas_obj, 'alpha', .5);
% surface_handles = isosurface(atlas_obj, 'alpha', .5, 'nomatchleftright');
%
% view(135, 30);
% lightRestoreSingle;
% lightFollowView;
%
% load(which('CIT168_MNI_subcortical_atlas_object.mat'));
%
% surface_handles = isosurface(atlas_obj, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]});
% surface_handles = isosurface(atlas_obj, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]}, 'nomatchleftright');
% p = addbrain('hires right');
% lightFollowView;

r = atlas2region(atlas_obj);

surface_handles = isosurface(r, varargin{:});


end % function


