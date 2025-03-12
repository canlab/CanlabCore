function [handles, r, subctx_fmri_data_obj] = render_cifti_on_brain(cifti_filename, varargin)
% render_cifti_on_brain renders an image from a cifti file onto surfaces and brain slices.
%
% :Usage:
% ::
%     render_cifti_on_brain(cifti_filename OR cifti_struct, [optional inputs])
%
% :Inputs:
%
%   **cifti_filename:**
%        String. Full path to the cifti file (e.g., 'transcriptomic_gradients.dscalar.nii').
%
% :Optional Inputs:
%
%   **'which_image':** [numeric scalar]
%        Index of the image (map) to render from the cifti file.
%        Default = 1.
%
%   **'color_map':** [matrix]
%        An n-by-3 colormap matrix used for rendering.
%        Default = colormap_tor([0 0 1], [1 1 0], [0.5 0.5 0.5], [0 0.5 1], [1 0.5 0]).
%
% :Outputs:
%
% **handles:**
%   graphics handles to axes and surfaces
%   use these to change colormaps and properties, etc.
%
% **r:**
% region object for subcortical volume structures
%
% **subctx_fmri_data_obj:**
% fmri_data object for subcortical volume structures
%
% :Dependencies:
% You need Wash U HCP CIFTI tools on your matlab path. 
% See https://github.com/Washington-University/HCPpipelines
%
% For cerebellar flatmap rendering, you need Diedrichsen's SUIT toolbox on
% your matlab path.
%
% :References:
%   Human Connectome Project cifti functions, Canlab core tools.
%
% :See also:
%   cifti_read, plot_surface_map, cifti_struct_2_region_obj
%   gifti, make_surface_figure,
%   colormap_tor, lightRestoreSingle, montage
%
% :Examples:
% ::
%     handles = render_cifti_on_brain('transcriptomic_gradients.dscalar.nii', 'which_image', 2, 'color_map', jet(256));
%
%     % Label the map
%     axes(handles.cortical_surface_axes(2));
%     title('                                         Transcriptomic gradient 2', 'Fontsize', 24)

% ..
%     Author and copyright information:
%
%     Copyright (C) 2025 Tor Wager
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

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------
%
% Some useful display helper functions
% n_cols = 80;
% sep_str = repmat('_', 1, n_cols);
% dashes = '----------------------------------------------';
% printstr = @(s) disp(s);
% printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

% Parse variable inputs using inputParser
ARGS = parse_inputs(varargin{:});

% Distribute parsed variables to workspace variables
% which_image = ARGS.which_image;
% color_map   = ARGS.color_map;

fn = fieldnames(ARGS);
for i = 1:length(fn)
    eval([fn{i}, ' = ARGS.(fn{i});']);
end

% Check for requirements
filename = which('cifti_struct_dense_extract_surface_data');
if isempty(filename), error('You need Wash U HCP CIFTI tools on your matlab path. See https://github.com/Washington-University/HCPpipelines'); end

% -------------------------------------------------------------------------
% Main function code: Render cifti data on brain surfaces and slices
% -------------------------------------------------------------------------

if isstruct(cifti_filename) && isfield(cifti_filename, 'metadata') && isfield(cifti_filename, 'cdata')
    % Assume this is a valid pre-loaded CIFTI structure
    cifti_struct = cifti_filename;

elseif ischar(cifti_filename)
    % Assume this is a CIFTI file name and load
    % Read the cifti file using HCP tools.
    cifti_struct = cifti_read(cifti_filename);

else
    error('cifti_filename is neither a valid CIFTI structure nor a filename')

end

% List available models
% "Models" are named structures with surface vertices or voxels
% --------------------------------------------------

model_names = listModelNames(cifti_struct);

if verbose
    disp('Models found in cifti file:');
    disp(model_names);
end

% % Get models array and number of models
% models = cifti_struct.diminfo{1}.models;

% Find surface models we will render specially 
wh_left = find(strcmp(model_names, 'CORTEX_LEFT'));
wh_right = find(strcmp(model_names, 'CORTEX_RIGHT'));


% Create a new figure for surface renderings.
% Set up surface axes
% --------------------------------------------------
create_figure('cifti surfaces');
set(gcf, 'Position', [20 groot().ScreenSize(4)*.25 groot().ScreenSize(3)*.6 groot().ScreenSize(4)*.6]);

xstart = 0.01;
cortexystart = 0.5; 
cortexwid = 0.95/4; 
spacegap = 0.01;

set(gca, 'Position', [xstart cortexystart cortexwid .95-cortexystart]);
ax(1) = gca;

ax(2) = axes('Position', [xstart+cortexwid+spacegap cortexystart cortexwid .95-cortexystart]);
ax(3) = axes('Position', [xstart+2*cortexwid+2*spacegap cortexystart cortexwid .95-cortexystart]);
ax(4) = axes('Position', [xstart+3*cortexwid+3*spacegap cortexystart cortexwid .95-cortexystart]);

% Load surfaces we will need
% --------------------------------------------------
l_surf = gifti('S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii'); % Load left cortical surface.
r_surf = gifti('S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii'); % Load right cortical surface.

% Render cortical colormap
% --------------------------------------------------

axes(ax(1))
han(1) = render_surface_patch(l_surf, color_map, [-90 0]);

axes(ax(2))
han(2) = render_surface_patch(l_surf, color_map, [-270 0]);

axes(ax(3))
han(3) = render_surface_patch(r_surf, color_map, [-90 0]);

axes(ax(4))
han(4) = render_surface_patch(r_surf, color_map, [-270 0]);

if wh_left
    leftdata = cifti_struct_dense_extract_surface_data(cifti_struct, 'CORTEX_LEFT', 1);
    leftdata = leftdata(:, which_image);

    set(han(1:2), 'FaceColor', 'interp', 'FaceVertexCData', leftdata);
end

if wh_right
    rightdata = cifti_struct_dense_extract_surface_data(cifti_struct, 'CORTEX_LEFT', 1);
    rightdata = rightdata(:, which_image);

    set(han(3:4), 'FaceColor', 'interp', 'FaceVertexCData', rightdata);
end

% get volume info as region object (r)
% --------------------------------------------------

r = cifti_struct_2_region_obj(cifti_struct, 'which_image', which_image);

% Show subcortical slices
% --------------------------------------------------
[o2, wh_montages] = make_slice_montage;

o2 = addblobs(o2, r, 'colormap', color_map, 'wh_montages', wh_montages);

% New figure with separate subcortical regions
% montage(r, 'colormap', color_map, 'regioncenters');

% Show subcortical surface
% --------------------------------------------------

sax_ystart = .28;

sax = axes('Position', [.22 sax_ystart .25 .25]);

subctx_han = addbrain('subcortex'); 
set(subctx_han, 'FaceAlpha', 1, 'FaceColor', [.5 .5 .5]); 
view(90, 1); 
lightRestoreSingle
surface(r, 'surface_handles', subctx_han, 'colormap', color_map, 'nolegend');
axis image

% surface() creates its own hybrid colormap to render gray and colored
% values separately, so we can't just use the one in color_map
mycm = get(sax(1), 'colormap');

% copy object to new view
sax(2) = axes('Position', [.45 sax_ystart .25 .25]);
subctx_han2 = copyobj(subctx_han, sax(2));
view(270, 0);
lightRestoreSingle
set(sax(2), 'colormap', mycm);

% re-set colormap of surfaces
% because montage has changed the default figure colormap
set(ax, 'colormap', color_map)
axis image
axis off

% Collect handles for output
% --------------------------------------------------
handles.cortical_surface_axes = ax;
handles.subcortical_surface_axes = sax;
handles.cortical_surfaces = han;
handles.subcortical_surfaces = [subctx_han subctx_han2'];
handles.fmri_display_o2_obj = o2;


% Show cerebellar flat map surface
% --------------------------------------------------
try
    subctx_fmri_data_obj = region2fmri_data(r);
catch
    disp('Error transforming region obj to fmri_data with region2fmri_data: fix me!')
    return
end

try

sax(3) = axes('Position', [.75 sax_ystart .25 .25]);
render_on_cerebellar_flatmap(subctx_fmri_data_obj, 'color_map', color_map);

set(sax(3), 'XLim', [-98 98], 'YLim', [-80 75]);
axis off

catch
    disp('Error rendering on cerebellum, render_on_cerebellar_flatmap. You may need SUIT and cb surface files on your matlab path')
    return
end


end % Main function

% transgrad = cifti_read('transcriptomic_gradients.dscalar.nii');
% r = cifti_struct_2_region_obj(cifti_struct, 'which_image', 3); % 3rd gradient
% subctx_fmri_data_obj = region2fmri_data(r);
% render_on_cerebellar_flatmap(subctx_fmri_data_obj, 'color_map', colormap('summer'))
% render_on_cerebellar_flatmap(subctx_fmri_data_obj)


% ------------------------------------------------------------------------
% Subfunction: parse_inputs
% ------------------------------------------------------------------------
function ARGS = parse_inputs(varargin)
% parse_inputs parses optional input arguments.
%
% :Usage:
% ::
%     ARGS = parse_inputs(optional_name_value_pairs)
%
% :Optional Inputs:
%
%   **'which_image':** [numeric scalar]
%        Index of the image (map) to render from the cifti file.
%        Default = 1.
%
%   **'color_map':** [matrix]
%        An n-by-3 colormap matrix used for rendering.
%        Default = colormap_tor([0 0 1], [1 1 0], [0.5 0.5 0.5], [0 0.5 1], [1 0.5 0]).
%
% ------------------------------------------------------------------------
p = inputParser;

addParameter(p, 'which_image', 1, @(x) isnumeric(x) && isscalar(x));

defaultColorMap = colormap_tor([0 0 1], [1 1 0], [0.5 0.5 0.5], [0 0.5 1], [1 0.5 0]);
addParameter(p, 'color_map', defaultColorMap, @(x) isnumeric(x) && size(x,2)==3);

addParameter(p, 'verbose', false, @(x) islogical(x) && isscalar(x));

parse(p, varargin{:});
ARGS = p.Results;

end

% ------------------------------------------------------------------------
% Subfunction: listModelNames
% ------------------------------------------------------------------------
function [model_names, model_types] = listModelNames(cifti_struct)
% listModelNames returns a cell array of model names from the cifti structure.

    models = cifti_struct.diminfo{1}.models;

    numModels = numel(models);
    [model_names, model_types] = deal(cell(numModels, 1));

    for i = 1:numModels
        model_names{i} = models{i}.struct;
        model_types{i} = models{i}.type;
    end

end % listModelNames

% ------------------------------------------------------------------------
% Subfunction: render_surface_patch
% ------------------------------------------------------------------------

function han = render_surface_patch(surface_struct, color_map, view_angle)

% han = patch('Vertices', surface_struct.vertices, 'Faces', surface_struct.faces, ...
%     'FaceColor', 'flat', 'EdgeColor', 'none');

han = patch('Faces',surface_struct.faces,'Vertices',surface_struct.vertices,'FaceColor', [.5 .5 .5], ...
'EdgeColor','none','SpecularStrength', .2, 'FaceAlpha', 1, 'SpecularExponent', 200);

material dull

colormap(color_map);

view(view_angle(1), view_angle(2));
axis off; axis image;
lightRestoreSingle;

end


% ------------------------------------------------------------------------
% Subfunction: make_slice_montage
% ------------------------------------------------------------------------


function [o2, wh_montages] = make_slice_montage

o2 = fmridisplay;

% set up custom montage
xstart = 0.2;
ystart = 0.05;
wid = 0.7/7;
spacegap = 0.0;

for i = 1:7
    ax(i) = axes('Position', [xstart+(i-1)*wid+(i-1)*spacegap ystart wid 2*wid]);
end

[o2, dat] = montage(o2, 'axial', 'slice_range', [-32 20], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', ax);

% set axes to zoom in
set(ax, 'XLim', [-50 50], 'YLim', [-100 55])

enlarge_axes(gcf, 1);
axh = axes('Position', [-0.06 .02 .29 .29]);  % [-0.02 0.15+shiftvals(i) .17 .17]);
axh(2) = axes('Position', [-0.02 .10 .29 .29]);

o2 = montage(o2, 'volume_data', dat, 'saggital', 'slice_range', [-2 2], 'spacing', 4, 'onerow', 'noverbose', 'existing_axes', axh);

wh_montages = [1 2];

brighten(.4)

end % make slice montage

