function [surface_handles, colors, patch_handles] = isosurface(atlas_obj, varargin)
% isosurface Create a series of surfaces in different colors, one per atlas region.
%
% Convert the atlas object to a region object and render each region as
% a 3-D isosurface. Supports custom colors and an option to disable
% automatic left/right color matching.
%
% :Usage:
% ::
%
%     [surface_handles, colors, patch_handles] = isosurface(atlas_obj, [optional arguments])
%
% :Inputs:
%
%   **atlas_obj:**
%        An atlas-class object.
%
% :Optional Inputs:
%
%   **Any optional inputs to imageCluster:**
%        e.g., 'alpha' followed by a value in [0,1] for surface
%        transparency.
%
%   **'colors':**
%        Followed by a single color in { } or a cell array of multiple
%        colors.
%
%   **'nomatchleftright':**
%        Do not match colors across hemispheres (left/right). The default
%        matches L/R and may override user-supplied colors.
%
% :Outputs:
%
%   **surface_handles:**
%        Vector of patch handles, one per region (legacy form; coerces
%        Patch objects to doubles).
%
%   **colors:**
%        Cell array of [r g b] color triplets used for each region.
%
%   **patch_handles:**
%        Cell array of patch object handles, one per region.
%
% :Examples:
% ::
%
%     atlasfile = which('Morel_thalamus_atlas_object.mat');
%     load(atlasfile)
%
%     surface_handles = isosurface(atlas_obj);
%     surface_handles = isosurface(atlas_obj, 'alpha', .5);
%     surface_handles = isosurface(atlas_obj, 'alpha', .5, 'nomatchleftright');
%
%     view(135, 30);
%     lightRestoreSingle;
%     lightFollowView;
%
%     load(which('CIT168_MNI_subcortical_atlas_object.mat'));
%
%     surface_handles = isosurface(atlas_obj, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]});
%     surface_handles = isosurface(atlas_obj, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]}, 'nomatchleftright');
%     p = addbrain('hires right');
%     lightFollowView;
%
% :See also:
%   - atlas2region
%   - region/isosurface
%   - imageCluster
%   - match_colors_left_right

r = atlas2region(atlas_obj);

[surface_handles, colors, patch_handles] = isosurface(r, varargin{:});


end % function


