function [surface_handles, pcl, ncl, o2] = surface(r, varargin)
% Render image data on brain surfaces; options for cutaways and canonical surfaces
% Surface method for region object - renders blobs on multiple types of 3-D
% surfaces
%
% :Usage:
% ::
%
%    [surface_handles, pcl, ncl, o2] = surface(r, [optional inputs])
%
%    - accepts combinations of existing surface handles and one or more
%    keywords accepted by addbrain.m to build surfaces
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2013 Tor Wager
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
% :Inputs:
%
%   **r:**
%        A region object
%
%   **'surface_handles':**
%        Followed by existing surface handles to render on
%
%   **'nocontroller':**
%        Render directly (legacy path) without registering the surfaces in an
%        fmridisplay object, so no interactive controller opens.
%
%   **addbrain.m keywords:**
%        Keywords for surfaces accepted by addbrain.m
%        A list is below:
%
%         % % Regions
%         % -----------------------------------------------------------------------
%         'vmpfc' 'nacc' 'BST' 'cau' 'caudate' 'put' 'GP' 'GPe' 'GPi' 'VeP' ...
%         'cm' 'md' 'stn' 'habenula' 'mammillary' 'hypothalamus','hy','hythal' ...
%         'midbrain' 'pag' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn' ...
%         'pbn' 'lc' 'rvm' 'rvm_old' 'nts' 'sc' 'ic' 'drn' 'mrn' ...
%         'thalamus' 'thal' 'LGN' 'lgn' 'MGN' 'mgn' 'VPthal', 'VPLthal', 'VPL', 'intralaminar_thal', ...
%         'medullary_raphe' 'spinal_trigeminal' 'nuc_ambiguus' 'dmnx_nts' 'ncs_B6_B8' 'nrp_B5' 'pbn' 'ncf' 'vep' 'PBP'
%         'amygdala' 'amygdala hires' 'hippocampus', 'hipp' 'hippocampus hires'
%
%         % Cortical Surfaces
%         % -----------------------------------------------------------------------
%         'left' 'hires left' 'surface left' 'hires surface left' ...
%         'right' 'hires right' 'surface right' 'hires surface right' ...
%         'transparent_surface' 'foursurfaces' 'foursurfaces_hcp' 'foursurfaces_hcp_inflated' 'flat left'  'flat right' ...
%
%         % Macro subcortical surfaces
%         % -----------------------------------------------------------------------
%         'CIT168' 'cerebellum','cblm' 'brainstem' 'suit brainstem'
%
%         % Cutaways
%         % -----------------------------------------------------------------------
%         'brainbottom' 'cutaway', 'left_cutaway' 'right_cutaway' ...
%         'left_insula_slab' 'right_insula_slab' 'accumbens_slab' 'coronal_slabs' 'coronal_slabs_4' 'coronal_slabs_5' ...
%
%         % Groups
%         % -----------------------------------------------------------------------
%         'bg', 'basal ganglia' 'midbrain_group' 'limbic' 'limbic hires' 'brainstem_group' 'thalamus_group'
%
%   **render_on_surface keywords:**
%   allowable = {'clim' 'color' 'colormap' 'colormapname' 'axis_handle' 'pos_colormap' 'neg_colormap'};
%   'color'         followed by [r g b] triplet for solid color
%   'clim'          Limits for image values mapped to most extreme colors
%   'colormap'      Followed by the name of a Matlab colormap to use (
%   'colormapname'
%   'axis_handle'
%   'pos_colormap'
%   'neg_colormap'
%
% :Outputs:
%
%   **all_surf_handles:**
%        surface patch handles
%
%   **pcl:**
%        region object with positive-only clusters
%
%   **ncl:**
%        region object with negative-only clusters
%
%   **o2:**
%        fmridisplay object holding the surfaces as registered views and the
%        region as a blob layer (empty with 'nocontroller'). The interactive
%        controller is bound to it, and you can keep working with it, e.g.
%        set_colormap(o2, 'color', [0 0 1]), removeblobs(o2), addblobs(o2, r2).
%
%
% :Examples:
% ::
%
%    % Render on combined cortical cutaway and subcortical surfaces
%    Use surface(r), with optional arguments taken by surface_cutaway:
%    poscm = colormap_tor([1 .3 0], [1 1 0]); % orange to yellow
%    [all_surf_handles, pcl, ncl] = surface(r, 'cutaway', 'ycut_mm', -30, 'pos_colormap', poscm, 'existingfig');
%
%    % Render on four cortical surfaces
%    [all_surf_handles2, pcl, ncl] = surface(r, 'foursurfaces', 'pos_colormap', poscm, 'neg_colormap', negcm);
%    [all_surf_handles2, pcl, ncl] = surface(r, 'foursurfaces', 'existingfig', 'color_upperboundpercentile', 95, 'color_lowerboundpercentile', 5, 'neg_colormap', colormap_tor([0 0 1], [.4 0 .7]));
%
%    % use mediation_brain_surface_figs and re-make colors [OLD METHOD! SLOW!]
%    all_surf_handles = mediation_brain_surface_figs([]);
%    surface(r, 'cutaway', 'surface_handles', all_surf_handles, 'color_upperboundpercentile', 95, 'color_lowerboundpercentile', 5, 'neg_colormap', colormap_tor([0 0 1], [.2 0 .5]));
%
%    % Make a region of interest surface and render on that
%    create_figure; p = addbrain('thalamus'); lightRestoreSingle;
%    [all_surf_handles2, pcl, ncl] = surface(region(t_age), 'color_upperboundpercentile', 95, 'color_lowerboundpercentile', 5, 'neg_colormap', colormap_tor([0 0 1], [.4 0 .7]), 'surface_handles', p);
%
% figure; hh = addbrain('foursurfaces');
% t.surface('surface_handles', hh);
%
% Note:
% To erase rendered blobs, use:
% sh = addbrain('eraseblobs', sh);
%
% :See also:*
%
% surface_cutaway, cluster_surf, mediation_brain_surface_figs

% ..
%    DEFAULTS AND INPUTS
% ..

existingfig = false;
surface_handles = [];
o2 = [];
use_display = true;   % register in an fmridisplay (enables the controller); 'nocontroller' -> legacy direct render

addbrain_allowable_args = {'vmpfc' 'nacc' 'BST' 'cau' 'caudate' 'put' 'GP' 'GPe' 'GPi' 'VeP' ...
    'cm' 'md' 'stn' 'habenula' 'mammillary' 'hypothalamus','hy','hythal' ...
    'midbrain' 'pag' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn' ...
    'pbn' 'lc' 'rvm' 'rvm_old' 'nts' 'sc' 'ic' 'drn' 'mrn' ...
    'thalamus' 'thal' 'LGN' 'lgn' 'MGN' 'mgn' 'VPthal', 'VPLthal', 'VPL', 'intralaminar_thal', ...
    'medullary_raphe' 'spinal_trigeminal' 'nuc_ambiguus' 'dmnx_nts' 'ncs_B6_B8' 'nrp_B5' 'pbn' 'ncf' 'vep' 'PBP' ...
    'left' 'hires left' 'surface left' 'hires surface left' ...
    'right' 'hires right' 'surface right' 'hires surface right' ...
    'transparent_surface' 'foursurfaces'  'foursurfaces_hcp' 'foursurfaces_hcp_inflated' 'flat left'  'flat right' 'foursurfaces_hcp' ...
    'brainbottom' 'cutaway', 'left_cutaway' 'right_cutaway' ...
    'left_insula_slab' 'right_insula_slab' 'accumbens_slab' 'coronal_slabs' 'coronal_slabs_4' 'coronal_slabs_5' ...
    'brainstem' 'suit brainstem' 'amygdala' 'amygdala hires' 'hippocampus', 'hipp' 'hippocampus hires' 'cerebellum','cblm' 'CIT168' ...
    'bg', 'basal ganglia' 'midbrain_group' 'limbic' 'limbic hires' 'brainstem_group' 'thalamus_group',...
    'inflated right','fsavg_right', 'inflated left', 'fsavg_left', 'hcp inflated left', 'hcp inflated right',...
    'hcp inflated', 'hcp sphere left', 'hcp sphere right', 'freesurfer sphere left', 'freesurfer sphere right',...
    'freesurfer inflated left', 'freesurfer inflated right', 'freesurfer white left', 'freesurfer white right',...
    'MNI152NLin2009cAsym white left','MNI152NLin2009cAsym white right','MNI152NLin2009cAsym midthickness left',...
    'MNI152NLin2009cAsym midthickness right','MNI152NLin2009cAsym pial left', 'MNI152NLin2009cAsym pial right',...
    'MNI152NLin6Asym white left','MNI152NLin6Asym white right','MNI152NLin6Asym midthickness left',...
    'MNI152NLin6Asym midthickness right','MNI152NLin6Asym pial left', 'MNI152NLin6Asym pial right',...
    'MNI152NLin6Asym sphere left','MNI152NLin6Asym sphere right','disableVis3d'};

render_on_surface_allowable_args = {'clim' 'cmaprange' 'color' 'colormap' 'colormapname' 'axis_handle' 'pos_colormap' 'neg_colormap' ...
    'summer' 'cool' 'hot' 'bone' 'copper' 'prism' 'hsv' 'winter' 'nolegend','sourcespace', 'targetsurface', 'MNI152NLin6Asym', 'fsLR_32k', 'MNI152NLin2009cAsym', 'fsaverage_164k', 'gray_buffer'};

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'surface_handles', surface_handles = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
            case {'existingfig', 'nofigure'}, existingfig = true; varargin{i} = []; 
                
            case 'noverbose'

            case 'nocontroller', use_display = false; varargin{i} = [];

            case addbrain_allowable_args
                % do nothing, handle later

            case render_on_surface_allowable_args
                % do nothing, will pass in to render_on_surface

            case {'unique', 'solid'}
                % colour-mode flags, translated to an indexed / solid colormap below

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% BUILD BRAIN SURFACE BASED ON ALLOWABLE INPUTS TO ADDBRAIN
% -------------------------------------------------------------------------

% Inflated / sphere / FreeSurfer surfaces are not in MNI space, so sampling
% the volume at their vertex coordinates paints the wrong place (usually
% almost nothing). If the caller didn't specify a registration, infer the
% target surface from the surface keyword, as region.montage does, and
% assume the CANlab default MNI152NLin6Asym (SPM) source space.
if ~any(strcmp(varargin, 'sourcespace')) && ~any(strcmp(varargin, 'targetsurface'))
    targetsurface = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if any(strcmp(varargin{i}, {'inflated left' 'inflated right' 'fsavg_left' 'fsavg_right'})) || ...
                    startsWith(varargin{i}, 'freesurfer ')
                targetsurface = 'fsaverage_164k';
            elseif any(strcmp(varargin{i}, {'hcp inflated left' 'hcp inflated right' 'hcp sphere left' 'hcp sphere right' 'foursurfaces_hcp_inflated'}))
                targetsurface = 'fsLR_32k';
            end
        end
    end
    if ~isempty(targetsurface)
        varargin = [varargin, {'sourcespace', 'MNI152NLin6Asym', 'targetsurface', targetsurface}];
    end
end

for i = 1:length(varargin)
    
    if ischar(varargin{i}) && any(strcmp(addbrain_allowable_args, varargin{i}))
        
        surface_handles = [surface_handles addbrain(varargin{i})];
        
    end
    
end

if isempty(surface_handles)
    
    surface_handles = addbrain('left_cutaway');
    
end


% END DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Legacy
if nargout > 1
    [pcl, ncl] = posneg_separate(r, 'Z');
end

% Run color change
% Convert to image vector first
% -------------------------------------------------------------------------
obj = region2imagevec(r);

if isempty(obj)
    disp('No data to display. Returning')
    return
end

if isempty(surface_handles)
    disp('No surface handles to display on. Returning')
    return
end

% Colour MODE: 'unique' (one solid colour per region, like region/montage),
% 'solid', or auto-selected via the shared canlab_color_mode when the caller gave
% no explicit colour spec. region2imagevec stores region indices in obj.dat, so
% 'unique' -> an indexed colormap (scn_standard_colors); render_on_surface uses
% nearest-neighbour automatically for an indexmap. Without this, an atlas rendered
% its integer region indices through a CONTINUOUS colormap (blended at borders).
explicit_color_keys = {'colormap', 'colormapname', 'color', 'pos_colormap', 'neg_colormap', 'indexmap'};
has_explicit_color = any(cellfun(@(a) (ischar(a) || isstring(a)) && any(strcmpi(a, explicit_color_keys)), varargin));
wu = find(strcmpi(varargin, 'unique'), 1);
ws = find(strcmpi(varargin, 'solid'), 1);
color_mode = '';
if ~isempty(wu),          color_mode = 'unique'; varargin(wu) = [];
elseif ~isempty(ws),      color_mode = 'solid';  varargin(ws) = [];
elseif ~has_explicit_color, color_mode = canlab_color_mode(obj);
end
switch color_mode
    case 'unique'
        nlab = max(1, round(max(double(obj.dat(:)))));
        varargin = [varargin, {'colormap', cell2mat(scn_standard_colors(nlab)'), 'indexmap', 'nolegend'}];
    case 'solid'
        varargin = [varargin, {'color', [1 0.4 0]}];
end

% Render through an fmridisplay object when possible: the surfaces become
% registered views and the region a managed blob layer, so the interactive
% controller opens and set_colormap / rethreshold / removeblobs work on it, as for
% montage(). The geometry above is unchanged. Falls back to a direct render.
if use_display
    try
        o2 = local_display_layer(r, surface_handles, varargin, addbrain_allowable_args);
    catch ME
        warning('region:surface:display', ...
            'Could not register surfaces in an fmridisplay (%s). Rendering directly.', ME.message);
        o2 = [];
    end
end
if isempty(o2)
    render_on_surface(obj, surface_handles, varargin{:});
end

% Hide the "..." interaction toolbar on the surface axes (recent MATLAB shows one
% in the corner of every axes). Sweep the figure(s) holding the surface patches.
try
    sh = surface_handles(ishandle(surface_handles));
    if ~isempty(sh)
        canlab_hide_axes_toolbar(ancestor(sh(1), 'figure'));
    end
catch
end

end % function


function o2 = local_display_layer(r, surface_handles, args, addbrain_allowable_args)
% Wrap already-built surface patches in an fmridisplay (one view per axes) and
% add region r as a layer. args are the remaining surface() options; addbrain
% keywords are geometry only, the rest are colour / registration options.

sh = surface_handles(ishandle(surface_handles));
if isempty(sh), error('no valid surface handles'); end

keep = true(size(args));
for i = 1:numel(args)
    if isempty(args{i}) || (ischar(args{i}) && any(strcmp(args{i}, [addbrain_allowable_args, {'noverbose'}])))
        keep(i) = false;
    end
end
args = args(keep);

% Translate render_on_surface colour conventions to addblobs' so the display
% path looks the same as the direct render:
%  - indexed ('unique') maps: render_on_surface takes 'colormap', M + an
%    'indexmap' flag; addblobs takes 'indexmap', M (nearest-neighbour throughout).
%  - no colour given: keep render_on_surface's hot/cool split rather than the
%    addblobs default palette.
wi = find(strcmp(args, 'indexmap'), 1);
wc = find(strcmp(args, 'colormap'), 1);
if ~isempty(wi) && ~isempty(wc) && isnumeric(args{wc + 1})
    M = args{wc + 1};
    args([wi, wc, wc + 1]) = [];
    args = [args, {'indexmap', M}];
end
color_keys = {'color', 'colormap', 'colormapname', 'pos_colormap', 'neg_colormap', 'splitcolor', 'maxcolor', 'mincolor', 'indexmap'};
if ~any(cellfun(@(a) ischar(a) && any(strcmp(a, color_keys)), args))
    args = [args, {'splitcolor', {[0 0 1], [0 .8 .8], [1 .4 .5], [1 1 0]}}];
end

o2 = fmridisplay;

% Group patches by parent axes, preserving creation order (as fmridisplay.surface
% does for addbrain's multi-panel layouts).
ax_list = gobjects(0);
grp = zeros(numel(sh), 1);
for i = 1:numel(sh)
    a = ancestor(sh(i), 'axes');
    idx = find(ax_list == a, 1);
    if isempty(idx), ax_list(end + 1) = a; idx = numel(ax_list); end %#ok<AGROW>
    grp(i) = idx;
end
for a = 1:numel(ax_list)
    these = sh(grp == a);
    o2.surface{end + 1} = struct('axis_handles', ax_list(a), 'direction', get(these(1), 'Tag'), ...
        'orientation', '', 'object_handle', these(:)');
end

o2 = addblobs(o2, r, args{:}, 'noverbose');

end

