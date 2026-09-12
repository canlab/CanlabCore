function o2 = canlab_montage_clusters(cl, varargin)
% Render clusters structures or region objects as colored layers on an fmridisplay slice montage
%
% :Usage:
% ::
%
%     o2 = canlab_montage_clusters(cl, [optional inputs])
%
% Modern replacement for the legacy montage_clusters and
% montage_clusters_medial functions. Accepts the old SCANlab "clusters"
% structure arrays as well as region objects, draws each set as a colored
% blob layer on a canlab_results_fmridisplay montage, and returns the
% fmridisplay object so more layers can be added with addblobs.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2026 Tor Wager
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
% ..
%
% :Inputs:
%
%   **cl:**
%        One of:
%        - a "clusters" structure array (legacy SCANlab format, e.g., from
%          mask2clusters or region2struct)
%        - a region object array
%        - a cell array of either, one cell per colored layer
%        Empty inputs (or empty cells) are allowed; the underlay is shown
%        with no blobs for that layer.
%
% :Optional Inputs:
%
%   **'colors', c:**
%        Cell array with one color per layer. Each entry can be a MATLAB
%        color name/char ('r', 'b', ...) or an [r g b] triplet. Colors are
%        recycled if fewer than the number of layers are given.
%        Default: {'r' 'b' 'g' 'c' 'm' 'y'}
%
%   **'colormap':**
%        Map the values stored in each region's .Z field (e.g., t or z
%        scores) to a color scale instead of using a solid color per layer.
%
%   **'overlay', img:**
%        Anatomical underlay image filename. Default: the standard
%        canlab_results_fmridisplay underlay.
%
%   **'montagetype', str:**
%        Any canlab_results_fmridisplay montage type, e.g. 'compact2'
%        [default], 'compact', 'full', 'multirow', 'coronal', 'sagittal'.
%        'sagittal' replaces the old montage_clusters_medial display.
%
%   **'o2', obj:**
%        An existing fmridisplay object. Layers are added to it and no new
%        figure is created.
%
%   **'figname', str:**
%        Name for the montage figure window.
%
%   **'addblobs_options', {...}:**
%        Cell array of extra options passed to addblobs for every layer,
%        e.g. {'outline', 'linewidth', 2} or {'trans', 'transvalue', .5}.
%
%   **'noverbose'** or **'verbose', false:**
%        Suppress text output.
%
% :Outputs:
%
%   **o2:**
%        fmridisplay object with one activation_maps entry per non-empty
%        layer.
%
% :Examples:
% ::
%
%    % Threshold a group t-map and display the regions
%    imgs = load_image_set('emotionreg', 'noverbose');
%    t = threshold(ttest(imgs), .001, 'unc', 'k', 10);
%    r = region(t);
%    o2 = canlab_montage_clusters(r, 'colors', {[1 .3 0]});
%
%    % Legacy clusters structure on sagittal slices
%    % (replaces montage_clusters_medial)
%    cl = region2struct(r);
%    o2 = canlab_montage_clusters(cl, 'colors', {'r'}, 'montagetype', 'sagittal');
%
%    % Two layers with different colors on a full montage
%    o2 = canlab_montage_clusters({r(1:3) r(4:end)}, 'colors', {'g' 'b'}, 'montagetype', 'full');
%
%    % Color blobs by their statistic values
%    o2 = canlab_montage_clusters(r, 'colormap', 'noverbose');
%
% :See also:
%   canlab_results_fmridisplay, region.montage, fmridisplay.addblobs,
%   cluster2region, region2struct

% ..
%    Programmers' notes:
%    2026-09  Created to replace montage_clusters and montage_clusters_medial
%             calls throughout the CANlab toolboxes.
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

% Parse special command keywords and remove them before inputParser

verbose = true;
verbose_idx = strcmpi(varargin, 'noverbose');
if any(verbose_idx)
    verbose = false;
    varargin(verbose_idx) = [];   % remove so inputParser doesn't see it
end

docolormap = false;
cmap_idx = strcmpi(varargin, 'colormap');
if any(cmap_idx)
    docolormap = true;
    varargin(cmap_idx) = [];      % remove so inputParser doesn't see it
end

% Use inputParser to parse key/value pairs
% First add obligatory/non-conditional keywords

p = inputParser;
p.addRequired('cl', @(x) isempty(x) || isstruct(x) || isa(x, 'region') || iscell(x));
p.addParameter('colors', {'r' 'b' 'g' 'c' 'm' 'y'}, @(x) iscell(x) || ischar(x) || isnumeric(x));
p.addParameter('overlay', '', @(x) isempty(x) || ischar(x) || isstring(x));
p.addParameter('montagetype', 'compact2', @(x) ischar(x) || isstring(x));
p.addParameter('o2', [], @(x) isempty(x) || isa(x, 'fmridisplay'));
p.addParameter('figname', '', @(x) ischar(x) || isstring(x));
p.addParameter('addblobs_options', {}, @iscell);

% Special key/value pairs that we have potentially set with optional keywords
p.addParameter('verbose', verbose, @(x) islogical(x) || isnumeric(x));

% process inputs and deal out to variables in workspace
p.parse(cl, varargin{:});

ARGS = p.Results;

colors = ARGS.colors;
overlay = char(ARGS.overlay);
montagetype = lower(char(ARGS.montagetype));
o2 = ARGS.o2;
figname = char(ARGS.figname);
addblobs_options = ARGS.addblobs_options;
verbose = logical(ARGS.verbose);

if ~iscell(colors), colors = {colors}; end
if ~iscell(cl), cl = {cl}; end

% canlab_results_fmridisplay spells this montage type 'saggital'
if strcmp(montagetype, 'sagittal'), montagetype = 'saggital'; end

% -------------------------------------------------------------------------
% Convert each layer to a region object
% -------------------------------------------------------------------------

nlayers = numel(cl);
layers = cell(1, nlayers);

for i = 1:nlayers
    layers{i} = to_region(cl{i});
end

nregions = sum(cellfun(@numel, layers));

if verbose
    fprintf('canlab_montage_clusters: %d layer(s), %d region(s), montage type ''%s''\n', nlayers, nregions, montagetype);
end

% -------------------------------------------------------------------------
% Create the underlay montage (unless an fmridisplay object was given)
% -------------------------------------------------------------------------

if isempty(o2)
    
    args = {montagetype, 'noblobs', 'nooutline'};
    if ~verbose, args{end + 1} = 'noverbose'; end
    if ~isempty(overlay), args = [args {'overlay', overlay}]; end
    
    o2 = canlab_results_fmridisplay(region(), args{:});
    
end

% -------------------------------------------------------------------------
% Add one blob layer per non-empty input
% -------------------------------------------------------------------------

ncolors = numel(colors);

for i = 1:nlayers
    
    if isempty(layers{i}), continue; end
    
    opts = addblobs_options;
    if ~verbose, opts{end + 1} = 'noverbose'; end %#ok<AGROW>
    
    if docolormap
        o2 = addblobs(o2, layers{i}, opts{:});
    else
        rgb = validatecolor(colors{1 + mod(i - 1, ncolors)});
        o2 = addblobs(o2, layers{i}, 'color', rgb, opts{:});
    end
    
end

if ~isempty(figname)
    fh = get_montage_figure(o2);
    if ~isempty(fh), set(fh, 'Name', figname); end
end

end % main function



% =========================================================================
% Subfunctions
% =========================================================================

function r = to_region(x)
% Convert a clusters structure (or region) to a region object; [] if empty

r = [];
if isempty(x), return; end

if isa(x, 'region')
    r = x;
elseif isstruct(x)
    r = cluster2region(x);
else
    error('canlab_montage_clusters:BadInput', ...
        'Each layer must be a clusters structure array or a region object.');
end

% Drop regions with no voxels
nvox = arrayfun(@(rr) size(rr.XYZmm, 2), r);
r = r(nvox > 0);

if isempty(r), r = []; end

end


function fh = get_montage_figure(o2)
% Figure handle that holds the first montage axes of an fmridisplay object

fh = [];
try
    fh = ancestor(o2.montage{1}.axis_handles(1), 'figure');
catch
end

if isempty(fh) && ~isempty(findall(0, 'Type', 'figure'))
    fh = gcf;
end

end
