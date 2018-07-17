function [all_surf_handles, pcl, ncl] = surface(r, varargin)
% Render image data on brain surfaces; options for cutaways and canonical surfaces
% Surface method for region object - renders blobs on multiple types of 3-D
% surfaces
%
% :Usage:
% ::
%
%    [all_surf_handles, pcl, ncl] = surface(r, ['cutaways', any optional inputs to surface_cutaway])
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
%   **'cutaway':**
%        String command for rendering cutaways instead of the default
%           - default is call to mediation_brain_surface_figs
%           - cutaways calls surface_cutaway
%           - all optional arguments are passed to surface_cutaway
%
%   **'rightsurface':**
%        String command for rendering a right frontal cortical
%                      view complementary to 'cutaways'
%
%   **'foursurfaces':**
%        Compact plots of four surfaces
%
%   **'surface_handles':**
%        Followed by surface handles to render on
%
% Other optional inputs to surface_cutaway
% e.g., 'pos_colormap', 'existingfig', 'mm_deep'
%
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
%    % use mediation_brain_surface_figs and re-make colors
%    all_surf_handles = mediation_brain_surface_figs([]);
%    surface(r, 'cutaway', 'surface_handles', all_surf_handles, 'color_upperboundpercentile', 95, 'color_lowerboundpercentile', 5, 'neg_colormap', colormap_tor([0 0 1], [.2 0 .5]));
%
%    % Make a region of interest surface and render on that
%    create_figure; p = addbrain('thalamus'); lightRestoreSingle;
%    [all_surf_handles2, pcl, ncl] = surface(region(t_age), 'color_upperboundpercentile', 95, 'color_lowerboundpercentile', 5, 'neg_colormap', colormap_tor([0 0 1], [.4 0 .7]), 'surface_handles', p);
%
%
% :See also:*
%
% surface_cutaway, cluster_surf, mediation_brain_surface_figs

% ..
%    DEFAULTS AND INPUTS
% ..

allowable_args = {'cutaway' 'rightsurface' 'existingfig', 'foursurfaces', 'surface_handles'}; % optional inputs
default_values = {0, 0, 0, 0, []};

% actions for inputs can be: 'assign_next_input' or 'flag_on'
actions = {'flag_on', 'flag_on', 'flag_on', 'flag_on', 'assign_next_input'};

% logical vector and indices of which inputs are text
textargs = cellfun(@ischar, varargin);
whtextargs = find(textargs);

for i = 1:length(allowable_args)
    
    % assign default
    % -------------------------------------------------------------------------
    
    eval([allowable_args{i} ' =  default_values{i};']);
    
    wh = strcmp(allowable_args{i}, varargin(textargs));
    
    if any(wh)
        % Optional argument has been entered
        % -------------------------------------------------------------------------
        
        wh = whtextargs(wh);
        if length(wh) > 1, warning(['input ' allowable_args{i} ' is duplicated.']); end
        
        switch actions{i}
            case 'assign_next_input'
                eval([allowable_args{i} ' = varargin{wh(1) + 1};']);
                
            case 'flag_on'
                eval([allowable_args{i} ' = 1;']);
                
            otherwise
                error(['Coding bug: Illegal action for argument ' allowable_args{i}])
        end
        
    end % argument is input
end

% END DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

[pcl, ncl] = posneg_separate(r, 'Z');

if ~isempty(surface_handles)
    %all_surf_handles = surface_handles;
    all_surf_handles = surface_cutaway('cl', r, varargin{:});
    
elseif cutaway
    all_surf_handles = surface_cutaway('cl', r, varargin{:});

elseif rightsurface
    % could easily be extended to any methods for addbrain.m, but may be
    % better to build your own and pass in custom handles.
    
    if ~existingfig, create_figure('right_surface'); end
    
    all_surf_handles = addbrain('hires');
    view(133, 10);
    set(all_surf_handles, 'FaceAlpha', .9);
    axis tight, axis image
    lightRestoreSingle; lighting gouraud
    
    all_surf_handles = surface_cutaway('cl', r, 'surface_handles', all_surf_handles, varargin{:});

elseif foursurfaces
    all_surf_handles = run_foursurfaces(r, existingfig, varargin{:});

else  %if ~(cutaway || rightsurface || foursurfaces)
    
    % default method
    all_surf_handles = mediation_brain_surface_figs({pcl}, {ncl});
    
end

% Or, could modify to use this as an option:
%cluster_surf(obj, varargin{:})

end






%% Subfunction

function all_surf_handles = run_foursurfaces(r, existingfig, varargin)

if ~existingfig, f1 = create_figure('foursurfaces'); else f1 = gcf; end

all_surf_handles = [];

nrows = 2;
ncols = 2;

% Right lateral
% ------------------------------------------------------------------------
figure(f1);
subplot(nrows, ncols , 1);
surfh = addbrain('hires right');
set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
view(90, 0)
lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

surfh2 = addbrain('brainstem');
surfh2 = [surfh2 addbrain('thalamus')];
set(surfh2, 'FaceColor', [.5 .5 .5], 'FaceAlpha', .8);
surfh = [surfh surfh2];

surfh = surface_cutaway('cl', r, 'surface_handles', surfh, varargin{:});

all_surf_handles = [all_surf_handles surfh];

% Right medial
% ------------------------------------------------------------------------

figure(f1);
axh = subplot(nrows, ncols , 4);
surfh2 = copyobj(surfh, axh);
view(270, 0);

all_surf_handles = [all_surf_handles surfh2'];

lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

% Left lateral
% ------------------------------------------------------------------------

figure(f1);
subplot(nrows, ncols , 2);
surfh = addbrain('hires left');
set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
view(270, 0)
lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

surfh2 = addbrain('brainstem');
surfh2 = [surfh2 addbrain('thalamus')];
set(surfh2, 'FaceColor', [.5 .5 .5], 'FaceAlpha', .8);
surfh = [surfh surfh2];

surfh = surface_cutaway('cl', r, 'surface_handles', surfh, varargin{:});

all_surf_handles = [all_surf_handles surfh];

% Left medial
% ------------------------------------------------------------------------

figure(f1);
axh = subplot(nrows, ncols , 3);
surfh2 = copyobj(surfh, axh);
view(90, 0);

all_surf_handles = [all_surf_handles surfh2'];

lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

end

