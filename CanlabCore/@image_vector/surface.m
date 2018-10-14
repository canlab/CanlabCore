function [all_surf_handles, pcl, ncl] = surface(obj, varargin)
% Render image data on brain surfaces; options for cutaways and canonical surfaces
%
% *Usage:*
%   [all_surf_handles, pcl, ncl] = surface(obj)
%   [all_surf_handles, pcl, ncl] = surface(r, ['cutaways', any optional inputs to surface_cutaway])
% 
% This function uses region.surface to create surface figures.
% See help region.surface for options.
%
% :Examples:
% ::
%
%    % Create an initial surface plot from an fmri_data object:
%    han = surface(regionmasks{2});  
%
%    % Now add a second region in green:
%    cluster_surf(region(regionmasks{2}), {[0 1 0]}, han, 5);
%
%    % Use optional arguments taken by surface_cutaway:
%    poscm = colormap_tor([1 .3 0], [1 1 0]); % orange to yellow
%    [all_surf_handles, pcl, ncl] = surface(t, 'cutaway', 'ycut_mm', -30, 'pos_colormap', poscm, 'existingfig');
%    [all_surf_handles2, pcl, ncl] = surface(t, 'foursurfaces', 'pos_colormap', poscm, 'neg_colormap', negcm);
%    [all_surf_handles2, pcl, ncl] = surface(t, 'foursurfaces', 'existingfig', 'color_upperboundpercentile', 95, 'color_lowerboundpercentile', 5, 'neg_colormap', colormap_tor([0 0 1], [.3 0 .5]));
%
%    % Use mediation_brain_surface_figs and re-make colors
%    all_surf_handles = mediation_brain_surface_figs([]);
%    surface(t2, 'cutaway', 'surface_handles', all_surf_handles, 'color_upperboundpercentile', 95, 'color_lowerboundpercentile', 5, 'neg_colormap', colormap_tor([0 0 1], [.2 0 .5]));
%

if size(obj.dat, 2) > 1
    obj = mean(obj);
end

if isa(obj, 'atlas')
    
    r = atlas2region(obj);
    
    
else
    
    r = region(obj);
    [all_surf_handles, pcl, ncl] = surface(r, varargin{:});

end


end

