function surface_handles = canlab_canonical_brain_surface_cutaways(method_keyword)
% Create one of a number of pre-set 3D brain views, and return surface handles for rendering
%
% :Usage:
% ::
%
% surface_handles = canlab_canonical_brain_surface_cutaways(method_keyword)
%
% - Uses fmri_data.isosurface method, which is a flexible way of creating surfaces
% - You can then render activation blobs on these surfaces using the
%   surface() object methods for CANlab objects or cluster_surf
% 
% - Uses a pre-set group-average anatomical image tuned for quality visual display
%   ('keuken_2014_enhanced_for_underlay.img') from Keuken et al. 7T atlas
%
%   Keuken et al. (2014). Quantifying inter-individual anatomical variability in the subcortex using 7T structural MRI
%   Forstmann et al. (2014). Forstmann, Birte U., Max C. Keuken, Andreas Schafer, Pierre-Louis Bazin, Anneke Alkemade, and Robert Turner. 2014. ?Multi-Modal Ultra-High Resolution Structural 7-Tesla MRI Data Repository.? Scientific Data 1 (December): 140050.
%
% - This function is called from addbrain, which has many options for
%   surface rendering.
%
% - fmri_data.isosurface() can be used to create many more custom surfaces
%   See this code for examples.
%
% :Inputs:
%
%   **method_keyword:** One of the following:
%     'left_cutaway'
%     'right_cutaway'
%     'left_insula_slab'
%     'right_insula_slab'
%     'accumbens_slab'
%
% :Outputs:
%   **A figure display with rendering**
%
%   **surface_handles:**
%        A vector of handles to surface objects (isosurfaces and isocaps)
%
% :See also:
%   fmri_data.isosurface, addbrain, cluster_cutaways, fmri_data.surface,
%   region.surface
%

% ..
%     Author and copyright information:
%
%     Copyright (C) <year>  <name of author>
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%

anat = fmri_data(which('keuken_2014_enhanced_for_underlay.img'), 'noverbose');

anat.source_notes = {};
anat.source_notes{1, 1} = 'Anatomical image: Adapted from the 7T high-resolution atlas of Keuken, Forstmann et al.';
anat.source_notes{2, 1} = 'Keuken et al. (2014). Quantifying inter-individual anatomical variability in the subcortex using 7T structural MRI. ';
anat.source_notes{3, 1} = 'Forstmann, Birte U., Max C. Keuken, Andreas Schafer, Pierre-Louis Bazin, Anneke Alkemade, and Robert Turner. 2014. ?Multi-Modal Ultra-High Resolution Structural 7-Tesla MRI Data Repository.? Scientific Data 1 (December): 140050.';

disp(char(anat.source_notes));

switch method_keyword
    
    case {'cutaway' 'cutaways' 'left_cutaway'}
        
        p = isosurface(anat, 'thresh', 140, 'nosmooth', 'ylim', [-Inf -30]);
        p2 = isosurface(anat, 'thresh', 140, 'nosmooth', 'xlim', [-Inf 0], 'YLim', [-30 Inf]);
        
        p3 = addbrain('limbic hires');
        alpha 1
        delete(p3(3)); p3(3) = [];
        set(p3, 'FaceAlpha', .6, 'FaceColor', [.5 .5 .5]);
        
        set_slab_view_properties(p, [91 0])
        
        surface_handles = [p p2 p3];
        
    case 'right_cutaway'
        
        p = isosurface(anat, 'thresh', 140, 'nosmooth', 'ylim', [-Inf -30]);
        p2 = isosurface(anat, 'thresh', 140, 'nosmooth', 'xlim', [0 Inf], 'yLim', [-30 Inf]);
        
        p3 = addbrain('limbic hires');
        alpha 1
        delete(p3(3)); p3(3) = [];
        set(p3, 'FaceAlpha', .6, 'FaceColor', [.5 .5 .5]);
        
        set_slab_view_properties(p, [-137 14])
        
        surface_handles = [p p2 p3];
        
    case 'left_insula_slab'
        
        p = isosurface(anat, 'thresh', 140, 'nosmooth', 'xlim', [15 44], 'yLim', [-45 Inf]);
        set_slab_view_properties(p, [91 0])
        surface_handles = p;
        
    case 'right_insula_slab'
        
        p = isosurface(anat, 'thresh', 140, 'nosmooth', 'xlim', [-15 -44], 'YLim', [-45 Inf]);
        set_slab_view_properties(p, [91 0])
        surface_handles = p;
        
    case 'accumbens_slab'
        
        p = isosurface(anat, 'thresh', 140, 'nosmooth', 'yLim', [-18 8]);
        set_slab_view_properties(p, [156 0])
        surface_handles = p;
        
    case 'coronal_slabs'
        
        surface_handles = render_coronal_slabs(anat, 6);
        
    case 'coronal_slabs_4'
        
        surface_handles = render_coronal_slabs(anat, 4);
        
    case 'coronal_slabs_5'
        
        surface_handles = render_coronal_slabs(anat, 5);
 
    otherwise
        error('illegal method keyword.');
        
end % switch

end % main function

function set_slab_view_properties(p, az_el)

set(p, 'FaceAlpha', 1);
colormap gray;
axis vis3d image
view(az_el(1), az_el(2));
lightRestoreSingle;

end


function p = render_coronal_slabs(anat, nslabs)

start_mm = -80;
end_mm = 50;
slab_thickness_mm = 20;

st = linspace(start_mm, end_mm - slab_thickness_mm, nslabs); % define starting points

axwid = 2*0.8 / (nslabs + 1);
xstart = .1;
ystart = .9 - axwid;
xpos = linspace(xstart, .65, nslabs);
ypos = linspace(ystart, .7 - axwid, nslabs);

p = {};

for i = 1:nslabs
    ax = axes('Position', [xpos(i) ypos(i) axwid axwid]);
    p{i} = isosurface(anat, 'thresh', 140, 'nosmooth', 'yLim', [st(i) st(i)+slab_thickness_mm]);
    set(ax, 'Color', 'none');
    set_slab_view_properties(p{i}, [135 12])
    axis off
    drawnow
end

p = cat(2, p{:});
end

