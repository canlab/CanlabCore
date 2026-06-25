function obj = surface(obj, varargin)
% Add one or more brain surfaces to an fmridisplay object (managed display).
%
% Each surface is registered as a view on the object, so blobs added with
% addblobs (and refresh / removeblobs / rethreshold) act on all surfaces and
% montages together. The actual surface geometry and default camera view come
% from addbrain, so ANY surface keyword addbrain understands works here.
%
% :Usage:
% ::
%
%     obj = surface(obj, 'direction', DIR, 'orientation', ORN, 'axes', POS)
%     obj = surface(obj, DIR)          % bare direction token (e.g. 'thalamus')
%
% :Inputs:
%
%   **obj:**
%        an fmridisplay object (handle). Initialize first, e.g. with montage.
%
% :Optional Inputs:
%
%   **'direction', DIR:**
%        Which surface to draw. DIR is passed through to addbrain, so the full
%        set of addbrain surface/region keywords is available (run `help
%        addbrain` for the complete, current list). DIR may also be given as a
%        bare token: surface(obj, 'thalamus') == surface(obj, 'direction',
%        'thalamus'). Common cortical choices:
%          'hires left' / 'hires right'              (default: 'hires right')
%          'inflated left' / 'inflated right'        (freesurfer/fsaverage)
%          'hcp inflated left' / 'hcp inflated right'
%          'freesurfer inflated left' / '... right'
%          'flat left' / 'flat right'
%        Composite / subcortical examples (all from addbrain):
%          'brainstem left' / 'brainstem right'      (brainstem + midbrain group)
%          'caudate left' / 'caudate right'          (basal ganglia group)
%          'thalamus', 'amygdala', 'cerebellum', 'limbic', 'bg', ...
%        Multi-surface keywords add a SET of views to this same object:
%          'foursurfaces', 'foursurfaces_hcp', 'foursurfaces_freesurfer'
%
%   **'orientation', ORN:**
%        'lateral' (default) or 'medial'. For an L/R cortical surface, 'medial'
%        mirrors the camera azimuth 180 degrees. Ignored for surfaces without a
%        meaningful medial view.
%
%   **'axes', POS:**
%        An existing axes handle, or a [left bottom width height] position
%        vector for a new axes. Default: current axes (gca).
%
%   Any remaining inputs are passed straight through to addbrain.
%
% :Outputs:
%
%   **obj:**
%        the fmridisplay object, with the new surface view(s) registered in
%        obj.surface and any existing blob layers rendered onto them.
%
% :Examples:
% ::
%
%     o2 = fmridisplay; o2 = montage(o2, 'axial');
%     o2 = surface(o2, 'direction', 'hires right', 'orientation', 'lateral');
%     o2 = surface(o2, 'thalamus');                 % bare addbrain keyword
%     o2 = surface(o2, 'foursurfaces_hcp');         % 2x2 set in its own figure
%
% :See also:
%   - addbrain (full surface/region keyword list), addblobs, montage, refresh,
%     removeblobs, render_on_surface, canlab_results_fmridisplay
%
% See help fmridisplay

if nargin == 0 || isempty(obj) || ~isa(obj, 'fmridisplay')
    obj = fmridisplay;
end
% initialize, if nothing passed in; but you would have to call overloaded
% method, fmridisplay.montage, to invoke this.


if isempty(obj.SPACE) || ~isstruct(obj.SPACE.V)
    error('fmridisplay is not initialized correctly. run obj = fmridisplay; first and then pass in to this method.')

end

% Multi-surface keywords (e.g. 'foursurfaces', 'foursurfaces_hcp') add a SET
% of surface views to THIS SAME object, laid out in the current figure. Each
% becomes its own registered view, so blobs/refresh act on all of them
% together. (fmridisplay is a handle class, so han is not replaced.)
multi = '';
for kw = {'foursurfaces_hcp', 'foursurfaces_freesurfer', 'foursurfaces'}
    if any(strcmp(varargin, kw{1})), multi = kw{1}; break; end
end
if ~isempty(multi)
    obj = add_multiple_surfaces(obj, multi, varargin{:});
    return
end

% Parse inputs. Three reserved keyword pairs ('direction', 'orientation',
% 'axes') are consumed here; the FIRST unrecognized char token is taken as a
% BARE direction (so surface(han, 'thalamus') works, not just surface(han,
% 'direction', 'thalamus')), and any remaining args pass straight through to
% addbrain. This is what lets ANY eligible addbrain surface keyword be used in
% the managed display (see `help addbrain` for the full surface/region list).
ax = gca;
dir = '';
orn = 'lateral';
extra = {};                                            % pass-through to addbrain
i = 1;
while i <= numel(varargin)
    a = varargin{i};
    if ischar(a)
        switch a
            case 'direction',   dir = varargin{i + 1}; i = i + 2; continue
            case 'orientation', orn = varargin{i + 1}; i = i + 2; continue
            case 'axes',        ax  = varargin{i + 1}; i = i + 2; continue
            otherwise
                if isempty(dir)
                    dir = a;                           % bare direction token
                else
                    extra{end + 1} = a; %#ok<AGROW>    % addbrain pass-through
                end
                i = i + 1; continue
        end
    else
        extra{end + 1} = a; %#ok<AGROW>                % non-char (e.g. addbrain flag value)
        i = i + 1;
    end
end
if isempty(dir), dir = 'hires right'; end

if isa(ax,'matlab.graphics.axis.Axes')
    axh = ax;
    axes(axh);
else
    axh = axes('Position', ax);
end

% Build the surface(s). addbrain owns the surface geometry AND the default
% (lateral) camera view + lighting for every keyword, including its many
% composites (brainstem/caudate groups, limbic, basal ganglia, ...), so this
% method is essentially a pass-through: ANY eligible addbrain surface works
% here automatically. We only (1) special-case 'bigbrain left/right', which map
% to a DIFFERENT addbrain surface than 'bigbrain left' would, and (2) mirror the
% azimuth 180 degrees for a medial view (medial = lateral azimuth + 180).
switch dir
    case 'bigbrain left'
        h = addbrain('bigbrain', extra{:});
        view(-137, 18); lightRestoreSingle;

    case 'bigbrain right'
        h = addbrain('bigbrain', extra{:});
        view(137, 18); lightRestoreSingle;

    otherwise
        try
            h = addbrain(dir, extra{:});
        catch err
            % addbrain rejects unknown keywords with a terse 'Unknown method.';
            % translate that into an informative message that names the offending
            % token and points at the full list (the rest of addbrain's errors,
            % e.g. a missing surface file, are passed through unchanged).
            if strcmp(err.message, 'Unknown method.')
                error('fmridisplay:surface:unknownDirection', ...
                    ['surface: ''%s'' is not a recognized addbrain surface/region keyword.\n' ...
                     'Run `help addbrain` for the full list. Common choices: ''hires left''/''hires right'', ' ...
                     '''inflated left''/''inflated right'', ''flat left''/''flat right''; cutaways are ' ...
                     '''left_cutaway''/''right_cutaway'' (note the order, not ''cutaway_left''); composites ' ...
                     'include ''insula surfaces'', ''brainstem left'', ''caudate right'', ''limbic''.'], dir);
            else
                rethrow(err);
            end
        end
        if strcmp(orn, 'medial')
            [az, el] = view(axh);
            view(axh, mod(az + 180, 360), el);
        end
end





if strcmp(h(1).FaceColor,'interp')
    set(h,'FaceAlpha', 1);
else
    set(h, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
end
%[
lightRestoreSingle; 
axis image;
axis off; 
lighting gouraud;
material dull;


% register info in fmridisplay object
% 'axis_handles', axh,
obj.surface{end + 1} = struct('axis_handles', axh, 'direction', dir, 'orientation', orn, 'object_handle', h);

% Pull-in: if blob layers already exist on this object, render them onto the
% surface we just added, so a surface added AFTER addblobs still shows the
% blobs and add/remove/refresh act on all views together (design notes 4.2).
new_surf_idx = numel(obj.surface);
for k = 1:numel(obj.activation_maps)
    obj = render_layer_surfaces(obj, k, new_surf_idx);
end

end % main function


function obj = add_multiple_surfaces(obj, multi, varargin)
% Add a canonical set of four surface views (L/R lateral + L/R medial) to the
% same fmridisplay object, arranged 2x2 in a dedicated figure. Each is added
% via a normal single-surface surface() call (handle semantics keep obj the
% same), which also triggers blob pull-in per surface.

% Ensure a dedicated figure for the 2x2 layout. If the current figure already
% has content (e.g. a montage), the four surfaces would otherwise be drawn over
% it (and appear to render nothing). Reuse an empty current figure (e.g. the
% user typed `figure;` first); otherwise open a new one.
curfig = get(groot, 'CurrentFigure');
if isempty(curfig) || ~isempty(findobj(curfig, 'Type', 'axes'))
    figure;
end
layout_fig = gcf;
set(layout_fig, 'Color', 'w');                       % white background behind the surfaces

rest = varargin;
rest(strcmp(rest, multi)) = [];                       % drop the multi keyword
for kw = {'direction', 'orientation', 'axes'}         % drop any single-view controls
    wh = find(strcmp(rest, kw{1}));
    for j = sort(wh, 'descend'), rest(j:j+1) = []; end
end

switch multi
    case 'foursurfaces_hcp'
        Lsurf = 'hcp inflated left';        Rsurf = 'hcp inflated right';
    case 'foursurfaces_freesurfer'
        Lsurf = 'freesurfer inflated left'; Rsurf = 'freesurfer inflated right';
    otherwise % 'foursurfaces'
        Lsurf = 'hires left';               Rsurf = 'hires right';
end

% {direction, orientation, axes-position} for a 2x2 layout:
% top row = lateral views, bottom row = medial views; left col = left hemi.
configs = { ...
    {Lsurf, 'lateral', [0.03 0.50 0.45 0.45]}, ...
    {Rsurf, 'lateral', [0.52 0.50 0.45 0.45]}, ...
    {Lsurf, 'medial',  [0.03 0.03 0.45 0.45]}, ...
    {Rsurf, 'medial',  [0.52 0.03 0.45 0.45]} };

for c = 1:numel(configs)
    cfg = configs{c};
    obj = surface(obj, 'direction', cfg{1}, 'orientation', cfg{2}, 'axes', cfg{3}, rest{:});
end

% Turn off the axis box/background for EVERY axes in the layout figure — the four
% surface axes, the hidden colorbar axes, and any stray default axes created
% during rendering (which otherwise shows a frame behind the surfaces). The
% surface patches are children of the axes and stay visible; colorbars are
% separate objects and are unaffected.
set(findobj(layout_fig, 'Type', 'axes'), 'Visible', 'off');

end

