function surface_handles = surface_cutaway(varargin)
% Make a specialized cutaway surface and add blobs if entered
%
% :Usage:
% ::
%
%    surface_handles = surface_cutaway(varargin)
%
% ..
%     Author and copyright information:
%     -------------------------------------------------------------------------
%     Copyright (C) 2013  Tor Wager
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
% :Optional Inputs:
%
%   With no inputs, this function creates a vector of surface handles and
%   returns them in surface_handles
%
%   **'cl':**
%        followed by clusters structure or region object with blob info
%
%   **'surface_handles':**
%        followed by vector of existing surface handles to plot blobs on
%
%   **'ycut_mm':**
%        followed by y-cutoff in mm for coronal section.
%          - if absent, shows entire medial surface
%
%   **'existingfig':**
%        Use existing figure axis; do not create new one
%          - if 'handles' are passed in, this is irrelevant
%
%   **'pos_colormap':**
%        followed by colormap for positive-going values
%          - n x 3 vector of rgb colors
%          - see colormap_tor or matlab's colormap
%
%   **'neg_colormap':**
%        followed by colormap for negative-going values
%
%   **'color_upperboundpercentile':**
%        followed by 1-100 percentile threshold; see Color values below
%
%   **'color_lowerboundpercentile':**
%        followed by 1-100 percentile threshold; see Color values below
%
% :Output:
%
%   **surface_handles:**
%        vector of surface handles for all surface objects
%
% :Color values:
%
% Creates color scale using color-mapped colors, thresholded based on
% percentile of the distribution of scores in .Z field.
%   - To change the upper bound for which percentile of Z-scores is mapped to
%     the highest color value, enter 'color_upperboundpercentile', [new percentile from 1-100]
%   - To change the lower bound value mapped to the lowest color value, 
%     enter 'color_lowerboundpercentile' followed by [new percentile from 1-100]
%
% :Examples:
% ::
%
%    % Create new surfaces:
%    % r is a region object made from a thresholded image (see region.m)
%    surface_handles = surface_cutaway('cl', r, 'ycut_mm', -30);
%
%    % Image blobs on existing handles:
%    surface_handles = surface_cutaway('cl', r, 'handles', surface_handles);
%
%    % Do it step by step:
%    p = surface_cutaway();
%    p = surface_cutaway('cl', r, 'surface_handles', p);
%
%    % Use custom colormaps (you can define 'pos_colormap' and 'neg_colormap'):
%    poscm = colormap_tor([.5 0 .5], [1 0 0]);  % purple to red
%    p = surface_cutaway('cl', r, 'surface_handles', p, 'pos_colormap', poscm);
%    p = surface_cutaway('cl', r, 'ycut_mm', -30, 'pos_colormap', poscm);
%
%    % use mediation_brain_surface_figs and re-make colors
%    all_surf_handles = mediation_brain_surface_figs([]);
%    surface(t2, 'cutaway', 'surface_handles', all_surf_handles, 'color_upperboundpercentile', 95, 'color_lowerboundpercentile', 5, 'neg_colormap', colormap_tor([0 0 1], [.2 0 .5]));


% ..
%    DEFAULTS AND INPUTS
% ..

poscm = colormap_tor([1 0 .5], [1 1 0], [.9 .6 .1]);  %reddish-purple to orange to yellow
negcm = colormap_tor([0 0 1], [0 1 1], [.5 0 1]);  % cyan to purple to dark blue
% Set color maps for + / - values

% optional inputs with default values
% -----------------------------------
% - allowable_args is a cell array of argument names
% - avoid spaces, special characters, and names of existing functions
% - variables will be assigned based on these names
%   i.e., if you use an arg named 'cl', a variable called cl will be
%   created in the workspace

allowable_args = {'cl', 'ycut_mm', 'pos_colormap', 'neg_colormap', ...
    'surface_handles', 'existingfig' 'color_upperboundpercentile', 'color_lowerboundpercentile' ...
    'mm_deep'};

default_values = {[], [], poscm, negcm, ...
    [], 0, 80, [], 4};

% define actions for each input
% -----------------------------------
% - cell array with one cell for each allowable argument
% - these have special meanings in the code below
% - allowable actions for inputs in the code below are: 'assign_next_input' or 'flag_on'

actions = {'assign_next_input', 'assign_next_input', 'assign_next_input', 'assign_next_input', ...
    'assign_next_input', 'flag_on', 'assign_next_input', 'assign_next_input', 'assign_next_input'};

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

if ~isempty(surface_handles)
    existingfig = true;
end

% Build brain surface
% -------------------------------------------------------------------------

if ~existingfig
    
    create_figure('brain surface');
    
end

if isempty(surface_handles)
    
    surface_handles = build_brain_surface(ycut_mm);
    
else
    
    whbad = ~ishandle(surface_handles);
    if any(whbad)
        warning('Some handles are not valid. Using those that are.');
        surface_handles(whbad) = [];
    end
    
end

if isempty(surface_handles)
    return
end

% quit if no cluster/region input

if isempty(cl), return, end

% add blobs
% -------------------------------------------------------------------------

% Set scale for colormap based on .Z field
% These are not necessarily Z-scores, but they are assumed to be intensity values for the voxels

% Color values:
% Creates color scale using color-mapped colors, thresholded based on
% percentile of the distribution of scores in .Z field.
% - To change the upper bound for which percentile of Z-scores is mapped to
% the highest color value, enter 'color_upperboundpercentile', [new percentile from 1-100]
% - To change the lower bound value mapped to the lowest color value, 
% enter 'color_lowerboundpercentile' followed by [new percentile from 1-100]
  
clZ = cat(2,cl.Z);

color_lowerboundvaluepos = 0;
color_lowerboundvalueneg = 0;

if ~isempty(color_lowerboundpercentile)
    color_lowerboundvaluepos = prctile(clZ(clZ > 0), color_lowerboundpercentile);
    color_lowerboundvalueneg = prctile(clZ(clZ < 0), 100-color_lowerboundpercentile);
end

% refZ = [min(clZ(clZ > 0)) max(clZ) min(clZ(clZ < 0)) min(clZ)];
% refZ = [0 prctile(clZ(clZ > 0), 80) prctile(clZ(clZ < 0), 20) 0];
refZ = [color_lowerboundvaluepos prctile(clZ(clZ > 0), color_upperboundpercentile) prctile(clZ(clZ < 0), 100-color_upperboundpercentile) color_lowerboundvalueneg];

% Add blobs to all surfaces
% - Map color scale to custom colormaps using 'heatmap' and refZ,
% - and also map intensity to transparency using 'colorscale'

% make bigger to avoid cross-hatching bug/problem
%specialp = false(size(surface_handles));
%specialp([end end-2]) = true;
cluster_surf(cl, mm_deep, 'heatmap', 'colormaps', pos_colormap, neg_colormap, surface_handles, refZ, 'colorscale');

%cluster_surf(cl, 4, 'heatmap', 'colormaps', pos_colormap, neg_colormap, surface_handles(~specialp), refZ, 'colorscale');

% create render_blobs_surface after
% [blobhan, cmaprange, mincolor, maxcolor] = render_blobs(currentmap, mymontage, SPACE, varargin)
% just take the input options, and create color map based on blobs.  and
% code separating positive and negative activation in some cases.  then
% just use one call to cluster_surf

end % main function


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function surface_handles = build_brain_surface(ycut_mm)

% BUILD BRAIN SURFACE

% handles
surface_handles = [];
p2 = [];
p3 = [];

surface_handles = addbrain('limbic hires');
delete(surface_handles(end)); surface_handles(end) = [];

set(surface_handles, 'FaceAlpha', .7)

surface_handles = [surface_handles addbrain('brainstem')];

set(gca, 'ZLim', [-85 82]);

set(surface_handles, 'FaceColor', [.5 .5 .5]);

hh = findobj(gcf, 'Tag', 'thalamus');
set(hh, 'FaceAlpha', .8);

hh = get(surface_handles, 'Tag');
% wh = strcmp(hh, 'hippocampus');
wh = strcmp(hh, 'HC');

set(surface_handles(wh), 'FaceAlpha', .6);

wh = strcmp(hh, 'caudate');
set(surface_handles(wh), 'FaceAlpha', .6);

% overlay = which('SPM8_colin27T1_seg.img');
% ovlname = 'SPM8_colin27T1_seg';
overlay = which('BigBrain_processed.nii'); %underlay...
ovlname = 'BigBrain_processed';

if ~isempty(ycut_mm)
    
    [D,Ds,hdr,p2,bestCoords] = tor_3d('whichcuts','y','coords',[0 ycut_mm 0], 'topmm', 90, 'filename', ovlname, 'intensity_threshold', 81);
    set(p2(1),'FaceColor',[.5 .5 .5]);
    
end

[D,Ds,hdr,p3,bestCoords] = tor_3d('whichcuts','x','coords',[-4 0 0], 'topmm', 90, 'filename', ovlname, 'intensity_threshold', 85, 'bottommm', -75);
set(p3(1),'FaceColor',[.5 .5 .5]);

    colormap(flipud(gray));
material dull;
axis off

view(133, 10);

lightRestoreSingle

lighting gouraud

surface_handles = [surface_handles p2 p3];

end





