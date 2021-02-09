function [surface_handles, pcl, ncl] = surface(r, varargin)
% Render image data on brain surfaces; options for cutaways and canonical surfaces
% Surface method for region object - renders blobs on multiple types of 3-D
% surfaces
%
% :Usage:
% ::
%
%    [surface_handles, pcl, ncl] = surface(r, [optional inputs])
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
%         'transparent_surface' 'foursurfaces' 'flat left'  'flat right' ...
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

addbrain_allowable_args = {'vmpfc' 'nacc' 'BST' 'cau' 'caudate' 'put' 'GP' 'GPe' 'GPi' 'VeP' ...
    'cm' 'md' 'stn' 'habenula' 'mammillary' 'hypothalamus','hy','hythal' ...
    'midbrain' 'pag' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn' ...
    'pbn' 'lc' 'rvm' 'rvm_old' 'nts' 'sc' 'ic' 'drn' 'mrn' ...
    'thalamus' 'thal' 'LGN' 'lgn' 'MGN' 'mgn' 'VPthal', 'VPLthal', 'VPL', 'intralaminar_thal', ...
    'medullary_raphe' 'spinal_trigeminal' 'nuc_ambiguus' 'dmnx_nts' 'ncs_B6_B8' 'nrp_B5' 'pbn' 'ncf' 'vep' 'PBP' ...
    'left' 'hires left' 'surface left' 'hires surface left' ...
    'right' 'hires right' 'surface right' 'hires surface right' ...
    'transparent_surface' 'foursurfaces' 'flat left'  'flat right' ...
    'brainbottom' 'cutaway', 'left_cutaway' 'right_cutaway' ...
    'left_insula_slab' 'right_insula_slab' 'accumbens_slab' 'coronal_slabs' 'coronal_slabs_4' 'coronal_slabs_5' ...
    'brainstem' 'suit brainstem' 'amygdala' 'amygdala hires' 'hippocampus', 'hipp' 'hippocampus hires' 'cerebellum','cblm' 'CIT168' ...
    'bg', 'basal ganglia' 'midbrain_group' 'limbic' 'limbic hires' 'brainstem_group' 'thalamus_group'};

render_on_surface_allowable_args = {'clim' 'color' 'colormap' 'colormapname' 'axis_handle' 'pos_colormap' 'neg_colormap' ...
    'summer' 'cool' 'hot' 'bone' 'copper' 'prism' 'hsv' 'winter'};

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'surface_handles', surface_handles = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
            case {'existingfig', 'nofigure'}, existingfig = true; varargin{i} = []; 
                
            case 'noverbose'
                
            case addbrain_allowable_args 
                % do nothing, handle later
                
            case render_on_surface_allowable_args
                % no nothing, will pass in to render_on_surface
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% BUILD BRAIN SURFACE BASED ON ALLOWABLE INPUTS TO ADDBRAIN
% -------------------------------------------------------------------------

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

render_on_surface(obj, surface_handles, varargin{:});

end % function