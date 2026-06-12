function varargout = canlab_orthviews(varargin)
% canlab_orthviews SPM-free three-panel orthogonal-slice viewer for CANlab objects.
%
% Drop-in replacement for the most-used spm_orthviews functionality, written
% without any SPM dependency. Renders a single registered figure with three
% tight panels (sagittal, coronal, axial) for an anatomical underlay, with
% interactive crosshair navigation by clicking, a right-click context menu
% for zooming, and a SPM-compatible API for adding/removing blob overlays.
%
% :Usage:
% ::
%
%     canlab_orthviews                                  % open default underlay
%     canlab_orthviews(image_obj)                       % fmri_data/statistic_image/atlas
%     canlab_orthviews(region_obj)                      % region(s)
%     canlab_orthviews(image_obj, 'overlay', file)      % custom underlay
%     canlab_orthviews(image_obj, 'posneg')             % red+/blue- solid blobs
%     canlab_orthviews(image_obj, 'color', {[0 1 0]})   % solid color blobs
%
%     xyz = canlab_orthviews('position')                % crosshair (mm) [1x3]
%     canlab_orthviews('Reposition', [10 -8 4])         % move crosshair
%     canlab_orthviews('zoom', 60)                      % set zoom radius mm
%     canlab_orthviews('zoom', Inf)                     % reset zoom
%     canlab_orthviews('xhairs', 'off')                 % hide crosshairs
%     canlab_orthviews('Reset')                         % close & clear state
%
%     canlab_orthviews('AddBlobs', [handle,] XYZ, Z, M)
%     canlab_orthviews('AddColouredBlobs', [handle,] XYZ, Z, M, col)
%     canlab_orthviews('AddBlobs', [handle,] obj [, opts...])  % obj = statistic_image / fmri_data / region / ...
%     canlab_orthviews('RemoveBlobs')                   % drop all blob layers
%
% Note: the SPM-style `handle` selects one of several linked orthviews
% windows. canlab_orthviews maintains a single linked figure, so the
% handle is accepted for source-compatibility but ignored — you can
% omit it entirely.
%
% :Inputs:
%
%   The first argument is dispatched by type:
%
%   **image_vector / fmri_data / statistic_image / atlas:**
%        The image is converted to one or more "blob" layers via region()
%        and rendered over the current underlay.
%
%   **region (array):**
%        Plotted directly as a blob layer.
%
%   **char/string:**
%        Treated as either a SPM-style command keyword (see the command
%        list above) or, if it is a path to a NIfTI file, as the underlay
%        image.
%
% :Optional Inputs:
%
%   **'overlay', file:**
%        Path / which()-resolvable filename of an anatomical underlay
%        image. Default: canlab_get_underlay_image (fmriprep20 template).
%
%   **'posneg':**
%        Plot positive values in warm colors and negative values in cool
%        colors as solid blobs.
%
%   **'color', cell:**
%        Cell array of one or more RGB triplets used as solid blob colors.
%
%   **'trans':**
%        Plot blobs with reduced alpha (semi-transparent).
%
%   **'doplot', [logical flag]:**
%        Create plots; default = true.
%
%   **'doverbose', [logical flag]:**
%        Verbose output; default = true.
%
% :Outputs:
%
%   When called as a command, the only side-effects-free command is
%   'position', which returns the current crosshair coordinate as a 1x3
%   row vector of millimetres. All other commands return nothing.
%
% :Examples:
% ::
%
%     % Display a thresholded t-map (matches the existing orthviews method)
%     dat   = load_image_set('emotionreg');
%     t     = ttest(dat);
%     t     = threshold(t, .005, 'unc', 'k', 10);
%     canlab_orthviews(t);
%
%     % Move the crosshair and read the position back
%     canlab_orthviews('Reposition', [-6 22 28]);
%     xyz = canlab_orthviews('position');
%
%     % Add a region object to the existing display
%     r = region(t);
%     canlab_orthviews(r);
%
% Example sequence of visualizations:
% obj = load_image_set('emotionreg');
% t = ttest(obj);
% canlab_orthviews(t, 'trans'); pause(2)
% canlab_orthviews('RemoveBlobs'); pause(2)
% t = threshold(t, .005, 'unc');
% canlab_orthviews('AddBlobs', t); pause(2)
% xyz = canlab_orthviews('position');
% canlab_orthviews('RemoveBlobs'); 
% canlab_orthviews('AddBlobs', t, 'color', {[0 1 0]})
% canlab_orthviews('Reposition', [57 -1 -23])
% canlab_orthviews('zoom', 20)
% canlab_orthviews('RemoveBlobs'); 
% canlab_orthviews('AddBlobs', t, 'unique')
%
% :See also:
%   - canlab_get_underlay_image
%   - cluster_orthviews
%   - canlab_results_fmridisplay
%   - @fmri_data/orthviews
%   - @region/orthviews
%
% ..
%    Copyright (C) 2026  Tor Wager / CANlab
%
%    Programmers' notes:
%    Image I/O uses MATLAB's built-in niftiread / niftiinfo (R2017b+).
%    Internal state is stored as appdata on the figure under the key
%    'canlab_orthviews_state'.
% ..

% -------------------------------------------------------------------------
% Top-level dispatch
% -------------------------------------------------------------------------

if nargin == 0
    % open the default underlay if no figure exists; otherwise just return
    [fh, state] = ensure_figure_with_underlay([]);
    if nargout > 0, varargout{1} = state.centre; end
    return
end

arg1 = varargin{1};

% Command-string dispatch: spm_orthviews-compatible commands
if (ischar(arg1) || (isstring(arg1) && isscalar(arg1))) ...
        && is_command_string(char(arg1))
    [varargout{1:nargout}] = dispatch_command(varargin{:});
    return
end

% Object/data dispatch: render the input on top of the underlay
[varargout{1:nargout}] = dispatch_object(varargin{:});

end % canlab_orthviews


% =========================================================================
% Command dispatch
% =========================================================================

function tf = is_command_string(s)
    cmds = {'reposition','pos','position','reset','xhairs','zoom', ...
            'image','addblobs','addcolouredblobs','addcoloredblobs', ...
            'removeblobs','rmblobs','rmcolouredblobs','rmcoloredblobs', ...
            'overlay','underlay','redraw','close','window','context', ...
            'legend', 'colorbar', 'nocolorbar','nolegend','no_colorbar','no_legend'};
    tf = ismember(lower(strtrim(s)), cmds);
end


function varargout = dispatch_command(varargin)
    cmd  = lower(strtrim(char(varargin{1})));
    args = varargin(2:end);
    varargout = {};

    switch cmd
        case {'reposition'}
            xyz = reshape(args{1}, 1, []);
            if numel(xyz) < 3, xyz(end+1:3) = 0; end
            [fh, state] = ensure_figure_with_underlay([]);
            state.centre = xyz(1:3);
            setappdata(fh, 'canlab_orthviews_state', state);
            redraw_all(fh);

        case {'pos','position'}
            fh = find_canlab_orthviews_figure();
            if isempty(fh)
                [fh, state] = ensure_figure_with_underlay([]);
            else
                state = getappdata(fh, 'canlab_orthviews_state');
            end
            varargout{1} = state.centre;

        case {'reset','close'}
            % Find and close every canlab_orthviews figure (defensively
            % handles the case where more than one window is open).
            f = findobj('Tag', 'canlab_orthviews');
            if ~isempty(f)
                is_fig = arrayfun(@(h) contains(class(h), 'Figure'), f);
                f = f(is_fig);
                if ~isempty(f), close(f); end
            end

        case 'xhairs'
            fh = find_canlab_orthviews_figure();
            if isempty(fh), return; end
            state = getappdata(fh, 'canlab_orthviews_state');
            if ~isempty(args) && ischar(args{1})
                state.show_xhairs = strcmpi(args{1}, 'on');
            else
                state.show_xhairs = ~state.show_xhairs;
            end
            setappdata(fh, 'canlab_orthviews_state', state);
            redraw_crosshairs(fh);

        case 'zoom'
            [fh, state] = ensure_figure_with_underlay([]);
            if isempty(args) || isempty(args{1})
                state.zoom = Inf;
            else
                state.zoom = double(args{1});
            end
            setappdata(fh, 'canlab_orthviews_state', state);
            redraw_all(fh);

        case {'overlay','underlay','image'}
            fname = args{1};
            [fh, ~] = ensure_figure_with_underlay(fname);
            % blobs already on the figure will be redrawn over the new underlay
            redraw_all(fh);

        case 'redraw'
            fh = find_canlab_orthviews_figure();
            if ~isempty(fh), redraw_all(fh); end

        case 'addblobs'
            % canlab_orthviews('AddBlobs', wh_handle, XYZ, Z, M)
            add_blob_layer(args, 'autocolor');

        case {'addcolouredblobs','addcoloredblobs'}
            % canlab_orthviews('AddColouredBlobs', wh_handle, XYZ, Z, M, col)
            add_blob_layer(args, 'solidcolor');

        case {'removeblobs','rmblobs','rmcolouredblobs','rmcoloredblobs'}
            fh = find_canlab_orthviews_figure();
            if isempty(fh), return; end
            state = getappdata(fh, 'canlab_orthviews_state');
            state.blobs = {};
            setappdata(fh, 'canlab_orthviews_state', state);
            redraw_all(fh);

        case 'window'
            % canlab_orthviews('window', [lo hi]) — set underlay window
            if ~isempty(args)
                fh = find_canlab_orthviews_figure();
                if isempty(fh), return; end
                state = getappdata(fh, 'canlab_orthviews_state');
                state.underlay.window = args{1}(:)';
                setappdata(fh, 'canlab_orthviews_state', state);
                redraw_all(fh);
            end

        case {'colorbar','legend'}
            set_colorbar_visibility(true);

        case {'nocolorbar','nolegend','no_colorbar','no_legend'}
            set_colorbar_visibility(false);

        case 'context'
            % no-op placeholder for parity with spm_orthviews('context_menu')
    end
end


function set_colorbar_visibility(showit)
    fh = find_canlab_orthviews_figure();
    if isempty(fh)
        [fh, ~] = ensure_figure_with_underlay([]);
    end
    state = getappdata(fh, 'canlab_orthviews_state');
    state.show_colorbar = logical(showit);
    setappdata(fh, 'canlab_orthviews_state', state);
    redraw_all(fh);
end


% =========================================================================
% Object dispatch
% =========================================================================

function varargout = dispatch_object(varargin)
    obj  = varargin{1};
    args = varargin(2:end);
    varargout = {};

    % Parse options that affect the underlay / coloring up front
    [overlay_file, opts, args] = parse_render_options(args);

    % Initialize / reuse figure with chosen underlay
    [fh, state] = ensure_figure_with_underlay(overlay_file);

    % Apply optional cosmetic overrides before rendering
    if ~isempty(opts.bgcolor)
        state.bgcolor = opts.bgcolor;
        setappdata(fh, 'canlab_orthviews_state', state);
    end
    if ~isempty(opts.show_xhairs)
        state.show_xhairs = logical(opts.show_xhairs);
        setappdata(fh, 'canlab_orthviews_state', state);
    end
    if ~isempty(opts.show_colorbar)
        state.show_colorbar = logical(opts.show_colorbar);
        setappdata(fh, 'canlab_orthviews_state', state);
    end

    % Drop existing blobs only when explicitly requested
    if opts.replace_blobs
        state.blobs = {};
        setappdata(fh, 'canlab_orthviews_state', state);
    end

    % Empty object input means "options only" — apply overrides then
    % redraw without adding blobs.
    if isempty(obj) || (isnumeric(obj) && isempty(obj))
        redraw_all(fh);
        if nargout > 0, varargout{1} = {}; end
        return
    end

    % Convert input object into one or more blob layers
    blob_layers = object_to_blobs(obj, opts);

    if isempty(blob_layers)
        redraw_all(fh);
        if nargout > 0, varargout{1} = {}; end
        return
    end

    state = getappdata(fh, 'canlab_orthviews_state');
    for k = 1:numel(blob_layers)
        state.blobs{end+1} = blob_layers{k};
    end
    setappdata(fh, 'canlab_orthviews_state', state);

    % Re-center on the largest blob's center of mass (matches SPM convention)
    state = recenter_on_first_blob(state);
    setappdata(fh, 'canlab_orthviews_state', state);

    redraw_all(fh);

    if nargout > 0, varargout{1} = blob_layers; end
end


function [overlay_file, opts, rest] = parse_render_options(args)
    overlay_file = '';
    opts = struct('posneg',false, 'colors',{{}}, ...
                  'transparency',0, ...        % 0 = opaque, 1 = fully transparent
                  'solid',false, 'replace_blobs',false, ...
                  'doverbose',true, ...
                  'bgcolor',[], ...            % [] = leave figure background unchanged
                  'show_xhairs',[], ...        % [] = leave crosshair state unchanged
                  'smooth_edges',true, ...     % anti-alias blob edges via trilinear
                  'smooth_sigma_mm',0, ...     % optional Gaussian pre-blur of blob volume
                  'unique',false, ...          % one solid-color layer per contiguous region
                  'show_colorbar',[]);         % [] = leave colorbar state unchanged

    rest = {};
    i = 1;
    while i <= numel(args)
        a = args{i};
        if ischar(a) || (isstring(a) && isscalar(a))
            key = lower(char(a));
            switch key
                case {'overlay','underlay'}
                    overlay_file = args{i+1};
                    i = i + 2; continue;
                case 'posneg'
                    opts.posneg = true;
                case 'color'
                    c = args{i+1};
                    if ~iscell(c), c = {c}; end
                    opts.colors = c;
                    i = i + 2; continue;
                case {'trans','transparent','transparency'}
                    % key-value pair: <0..1>, where 0 is opaque and 1 is
                    % fully transparent. Bare keyword keeps the legacy
                    % "moderately transparent" default (~0.4) so old
                    % 'trans' calls continue to look the same.
                    if i + 1 <= numel(args) && isnumeric(args{i+1}) ...
                            && isscalar(args{i+1}) ...
                            && args{i+1} >= 0 && args{i+1} <= 1
                        opts.transparency = double(args{i+1});
                        i = i + 2; continue;
                    else
                        opts.transparency = 0.4;
                    end
                case 'unique'
                    opts.unique = true;
                case {'nocolorbar','nolegend','no_colorbar','no_legend'}
                    opts.show_colorbar = false;
                case {'colorbar','legend'}
                    opts.show_colorbar = true;
                case 'solid'
                    opts.solid = true;
                case {'replaceblobs','clear','new'}
                    opts.replace_blobs = true;
                case 'noverbose'
                    opts.doverbose = false;
                case {'backgroundcolor','bgcolor'}
                    opts.bgcolor = resolve_color(args{i+1});
                    i = i + 2; continue;
                case {'black','blackbg'}
                    opts.bgcolor = [0 0 0];
                case {'white','whitebg'}
                    opts.bgcolor = [1 1 1];
                case 'crosshairs'
                    % key-value pair 'crosshairs', <true/false>; if the next
                    % argument isn't a boolean/numeric scalar, treat the
                    % keyword alone as "show crosshairs"
                    consumed = false;
                    if i + 1 <= numel(args)
                        n = args{i+1};
                        if islogical(n) || (isnumeric(n) && isscalar(n))
                            opts.show_xhairs = logical(n);
                            consumed = true;
                        end
                    end
                    if ~consumed, opts.show_xhairs = true; end
                    if consumed, i = i + 2; continue; end
                case 'nocrosshairs'
                    opts.show_xhairs = false;
                case {'smooth_edges','smoothedges'}
                    opts.smooth_edges = true;
                case {'no_smooth_edges','nosmoothedges','no_smooth','nosmooth'}
                    opts.smooth_edges = false;
                case 'smooth'
                    % 'smooth', sigma_mm — Gaussian pre-blur of blob volumes
                    if i + 1 <= numel(args) && isnumeric(args{i+1}) ...
                            && isscalar(args{i+1})
                        opts.smooth_sigma_mm = max(0, double(args{i+1}));
                        i = i + 2; continue;
                    else
                        % bare 'smooth' keyword -> modest default 2 mm
                        opts.smooth_sigma_mm = 2;
                    end
                otherwise
                    rest{end+1} = a; %#ok<AGROW>
            end
        else
            rest{end+1} = a; %#ok<AGROW>
        end
        i = i + 1;
    end
end


function rgb = resolve_color(c)
    % Accept a 3-element RGB (0-1 or 0-255), a MATLAB color letter, or a
    % common color name.
    if isnumeric(c) && numel(c) >= 3
        rgb = double(c(1:3));
        if any(rgb > 1), rgb = rgb / 255; end
        return
    end
    if ischar(c) || (isstring(c) && isscalar(c))
        switch lower(char(c))
            case {'white','w'},  rgb = [1 1 1];
            case {'black','k'},  rgb = [0 0 0];
            case {'red','r'},    rgb = [1 0 0];
            case {'green','g'},  rgb = [0 0.85 0];
            case {'blue','b'},   rgb = [0 0 1];
            case {'yellow','y'}, rgb = [1 1 0];
            case {'cyan','c'},   rgb = [0 1 1];
            case {'magenta','m'},rgb = [1 0 1];
            case {'gray','grey'},rgb = [0.5 0.5 0.5];
            otherwise,           rgb = [1 1 1];
        end
    else
        rgb = [1 1 1];
    end
end


function c = fg_for(bg)
    % Pick a foreground (text/axis label) color with reasonable contrast
    % against the given background.
    if isempty(bg), bg = [1 1 1]; end
    if mean(bg(1:3)) > 0.5, c = [0 0 0]; else, c = [1 1 1]; end
end


function layers = object_to_blobs(obj, opts)
    layers = {};

    if isa(obj, 'region')
        layers = region_to_layers(obj, opts);

    elseif isa(obj, 'image_vector') || isa(obj, 'fmri_data') ...
            || isa(obj, 'statistic_image') || isa(obj, 'atlas')
        layers = imgobj_to_layers(obj, opts);

    elseif ischar(obj) || isstring(obj)
        try
            d = fmri_data(char(obj));
            layers = imgobj_to_layers(d, opts);
        catch
            warning('canlab_orthviews:input', ...
                'Could not interpret input "%s" as a NIfTI/Analyze image.', char(obj));
        end

    else
        warning('canlab_orthviews:input', ...
            'Unrecognized input class: %s', class(obj));
    end
end


function layers = imgobj_to_layers(obj, opts)
    % Build blob layers from an image_vector / fmri_data / statistic_image /
    % atlas. Uses reconstruct_image so statistic_image .sig thresholds are
    % honored automatically (reconstruct_image multiplies .dat by .sig).
    obj = replace_empty(obj);

    % atlas: integer-label colormap, one layer per label
    if isa(obj, 'atlas')
        if isa(obj, 'statistic_image')
            V = reconstruct_image(obj);
        else
            V = iimg_reconstruct_vols(double(obj.dat), obj.volInfo);
        end
        if ndims(V) > 3, V = V(:, :, :, 1); end
        layers = {make_blob_layer(V, obj.volInfo.mat, 'unique', [], opts)};
        return
    end

    % 'unique': one solid-color layer per contiguous cluster (matches the
    % legacy image_vector.orthviews 'unique' behavior).
    if opts.unique
        layers = imgobj_unique_layers(obj, opts);
        return
    end

    if isa(obj, 'statistic_image')
        V = reconstruct_image(obj);
    else
        V = iimg_reconstruct_vols(double(obj.dat), obj.volInfo);
    end
    if ndims(V) > 3
        V = V(:, :, :, 1);
    end
    M = obj.volInfo.mat;

    if opts.posneg
        layers = split_posneg(V, M, opts);
        return
    end

    if ~isempty(opts.colors)
        layers = {make_blob_layer(V, M, 'solid', opts.colors{1}, opts)};
        return
    end

    layers = {make_blob_layer(V, M, 'autocolor', [], opts)};
end


function layers = imgobj_unique_layers(obj, opts)
    % Split an image_vector-derived object into contiguous regions and
    % build one solid-color blob layer per region. Colors come from
    % scn_standard_colors, matching @image_vector/orthviews 'unique'.
    obj = replace_empty(obj);
    try
        r = region(obj, 'contiguous_regions');
    catch ME
        warning('canlab_orthviews:unique', ...
            'region(obj, ''contiguous_regions'') failed: %s', ME.message);
        layers = {};
        return
    end
    layers = region_unique_layers(r, obj.volInfo.mat, obj.volInfo.dim, opts);
end


function layers = region_unique_layers(r, M, dim, opts)
    % One solid-color blob layer per region in the array.
    layers = {};
    if isempty(r), return; end
    if isempty(M), M = eye(4); end
    if nargin < 4 || isempty(dim)
        XYZall = [];
        for k = 1:numel(r)
            if ~isempty(r(k).XYZ)
                XYZall = [XYZall r(k).XYZ];                 %#ok<AGROW>
            end
        end
        if isempty(XYZall), return; end
        dim = max(round(XYZall), [], 2)';
        dim = max(dim, [1 1 1]);
    end

    try
        colors = scn_standard_colors(numel(r));
    catch
        colors = arrayfun(@(k) [rand rand rand], 1:numel(r), ...
            'UniformOutput', false);
    end

    for k = 1:numel(r)
        if isempty(r(k).XYZ), continue; end
        XYZ = round(r(k).XYZ);
        XYZ(XYZ < 1) = 1;
        keep = XYZ(1, :) <= dim(1) & XYZ(2, :) <= dim(2) & XYZ(3, :) <= dim(3);
        XYZ = XYZ(:, keep);
        if isempty(XYZ), continue; end

        zk = r(k).Z;
        if isempty(zk) || numel(zk) < size(r(k).XYZ, 2)
            zk = ones(1, size(r(k).XYZ, 2));
        end
        zk = zk(keep);
        if isempty(zk), zk = ones(1, size(XYZ, 2)); end

        V = zeros(dim, 'single');
        V(sub2ind(dim, XYZ(1, :), XYZ(2, :), XYZ(3, :))) = single(zk);

        layers{end+1} = make_blob_layer(V, M, 'solid', colors{k}, opts); %#ok<AGROW>
    end
end


function layers = region_to_layers(r, opts)
    if isempty(r), layers = {}; return; end
    M = r(1).M;
    if isempty(M), M = eye(4); end

    % 'unique': one solid-color layer per region in the array.
    if opts.unique
        layers = region_unique_layers(r, M, [], opts);
        return
    end

    % Rasterize all region voxel indices into a dense volume.
    XYZ = []; Z = [];
    for k = 1:numel(r)
        if isempty(r(k).XYZ), continue; end
        XYZ = [XYZ r(k).XYZ];                       %#ok<AGROW>
        zk = r(k).Z;
        if isempty(zk) || numel(zk) < size(r(k).XYZ, 2)
            zk = ones(1, size(r(k).XYZ, 2));
        end
        Z = [Z reshape(zk, 1, [])];                 %#ok<AGROW>
    end
    if isempty(XYZ), layers = {}; return; end

    XYZ = round(XYZ);
    XYZ(XYZ < 1) = 1;
    dim = max(XYZ, [], 2)';
    dim = max(dim, [1 1 1]);

    V = zeros(dim, 'single');
    linidx = sub2ind(dim, XYZ(1, :), XYZ(2, :), XYZ(3, :));
    V(linidx) = single(Z);

    if opts.posneg
        layers = split_posneg(V, M, opts);
        return
    end

    if ~isempty(opts.colors)
        layers = {make_blob_layer(V, M, 'solid', opts.colors{1}, opts)};
        return
    end

    layers = {make_blob_layer(V, M, 'autocolor', [], opts)};
end


function layers = split_posneg(V, M, opts)
    Vp = V; Vp(~(Vp > 0 & isfinite(Vp))) = 0;
    Vn = V; Vn(~(Vn < 0 & isfinite(Vn))) = 0;
    layers = {};
    if any(Vp(:) ~= 0)
        layers{end+1} = make_blob_layer(Vp, M, 'solid', [1 0.55 0], opts);
    end
    if any(Vn(:) ~= 0)
        layers{end+1} = make_blob_layer(Vn, M, 'solid', [0.2 0.4 1], opts);
    end
end


function L = make_blob_layer(V, M, style, color, opts)
    % Layer struct backed by a 3-D blob volume + voxel-to-mm affine. The
    % XYZmm/Z fields are derived for back-compat (recenter, tests) but the
    % renderer reads .vol so each underlying voxel covers all panel pixels
    % within its mm footprint (no more single-pixel dots).
    %
    % Smoothing knobs (set via opts in parse_render_options):
    %   smooth_sigma_mm > 0 — Gaussian pre-blur of the blob volume in mm
    %                         space. Expands the blob footprint smoothly.
    %   smooth_edges = true — trilinear sampling + alpha-from-fractional-
    %                         mask at render time, so blob edges fade
    %                         smoothly into the background.

    V = single(V);
    mask = single(V ~= 0 & isfinite(V));   % keep the original support

    % Optional Gaussian pre-blur of BOTH the values and the binary mask.
    % Smoothing them independently is important: smooth3(V) keeps a long
    % Gaussian tail (because the Gaussian's support is everywhere the
    % kernel touches), so deriving the mask from the smoothed values
    % would give a "dilated binary" alpha. By smoothing the mask itself
    % we get a true Gaussian alpha falloff, which is what makes the
    % rendered blob look smooth instead of swelling into big blocks.
    % (Skipped for atlas/unique style — would interpolate between
    % integer labels and destroy them.)
    if opts.smooth_sigma_mm > 0 && ~strcmp(style, 'unique')
        vox_mm = abs(det(M(1:3, 1:3)))^(1/3);
        if isfinite(vox_mm) && vox_mm > 0
            sigma_v = opts.smooth_sigma_mm / vox_mm;
            if sigma_v >= 0.5
                k = 2 * ceil(2 * sigma_v) + 1;     % cover ±2σ
                V    = single(smooth3(double(V),    'gaussian', k, sigma_v));
                mask = single(smooth3(double(mask), 'gaussian', k, sigma_v));
            end
        end
    end

    L = struct( ...
        'vol',             V, ...
        'mask',            mask, ...
        'M',               M, ...
        'dim',             size(V), ...
        'style',           style, ...
        'color',           color, ...
        'alpha',           1, ...
        'clim',            [0 1], ...
        'smooth_edges',    opts.smooth_edges, ...
        'smooth_sigma_mm', opts.smooth_sigma_mm);

    L.alpha = max(0, min(1, 1 - opts.transparency));

    nz = V(V ~= 0 & isfinite(V));
    if ~isempty(nz)
        L.clim = double([min(abs(nz)) max(abs(nz))]);
    end

    [L.XYZmm, L.Z] = volume_to_xyzmm(V, M);
end


function [XYZmm, Z] = volume_to_xyzmm(V, M)
    mask = V ~= 0 & isfinite(V);
    if ~any(mask(:))
        XYZmm = zeros(3, 0); Z = zeros(1, 0); return
    end
    [I, J, K] = ind2sub(size(V), find(mask));
    P = [I(:)'; J(:)'; K(:)'; ones(1, numel(I))];
    mm = M * P;
    XYZmm = mm(1:3, :);
    Z = double(V(mask))';
end


function add_blob_layer(args, mode)
    % Accepted call shapes:
    %   canlab_orthviews('AddBlobs', h, XYZ, Z, M)              % SPM-style
    %   canlab_orthviews('AddBlobs',    XYZ, Z, M)              % handle omitted
    %   canlab_orthviews('AddColouredBlobs', h, XYZ, Z, M, col) % SPM-style
    %   canlab_orthviews('AddColouredBlobs',    XYZ, Z, M, col) % handle omitted
    %   canlab_orthviews('AddBlobs', [h,] obj [, opts...])      % object input
    %
    % The numeric handle `h` exists in SPM to select one of several linked
    % orthviews windows. canlab_orthviews maintains a single linked
    % figure, so the handle is accepted for compatibility but ignored.
    % `obj` may be statistic_image / fmri_data / image_vector / atlas /
    % region (array); options like 'unique', 'transparent', 'color',
    % 'smooth', 'posneg' all work the same as in the direct object call.

    if isempty(args)
        warning('canlab_orthviews:addblobs', ...
            'AddBlobs requires at least one argument.');
        return
    end

    % Strip the optional SPM handle (a lone scalar number).
    if isnumeric(args{1}) && isscalar(args{1})
        args = args(2:end);
    end
    if isempty(args)
        warning('canlab_orthviews:addblobs', ...
            'AddBlobs needs a data argument.');
        return
    end

    first = args{1};

    % --- Object-style call ------------------------------------------------
    if isa(first, 'image_vector')   || isa(first, 'fmri_data') || ...
       isa(first, 'statistic_image') || isa(first, 'atlas')    || ...
       isa(first, 'region')
        % Forward to the standard dispatcher; all object options
        % ('unique', 'color', 'transparent', 'smooth', 'posneg', ...) are
        % parsed there. AddColouredBlobs with a missing color is a no-op
        % in this branch — the user passes color via the 'color' option.
        dispatch_object(first, args{2:end});
        return
    end

    % --- Legacy SPM-compatible XYZ / Z / M call --------------------------
    if numel(args) < 3
        warning('canlab_orthviews:addblobs', ...
            'AddBlobs(XYZ, Z, M) needs XYZ, Z, and M.');
        return
    end
    XYZ = args{1};
    Z   = args{2};
    M   = args{3};
    if size(M, 1) ~= 4 || size(M, 2) ~= 4
        warning('canlab_orthviews:addblobs', 'M must be 4x4.');
        return
    end

    XYZ = round(XYZ);
    XYZ(XYZ < 1) = 1;
    dim = max(XYZ, [], 2)';
    dim = max(dim, [1 1 1]);
    V = zeros(dim, 'single');
    linidx = sub2ind(dim, XYZ(1, :), XYZ(2, :), XYZ(3, :));
    V(linidx) = single(reshape(Z, 1, []));

    opts = struct('posneg',false, 'colors',{{}}, ...
                  'transparency',0, ...
                  'solid',false, 'replace_blobs',false, 'doverbose',true, ...
                  'smooth_edges', true, 'smooth_sigma_mm', 0, ...
                  'unique', false);

    if strcmp(mode, 'solidcolor') && numel(args) >= 4
        L = make_blob_layer(V, M, 'solid', args{4}, opts);
    else
        L = make_blob_layer(V, M, 'autocolor', [], opts);
    end

    [fh, state] = ensure_figure_with_underlay([]);
    state.blobs{end+1} = L;
    setappdata(fh, 'canlab_orthviews_state', state);
    redraw_all(fh);
end


% =========================================================================
% Figure setup and underlay loading
% =========================================================================

function [fh, state] = ensure_figure_with_underlay(overlay_file)
    fh = find_canlab_orthviews_figure();
    if isempty(fh)
        fh = build_figure();
    else
        % Refresh interaction wiring on every entry. This replaces any
        % stale callback handles that an older version of this file may
        % have attached, which is what causes "Unable to find function
        % @(src,evt)on_axis_click(...)" on reused figures.
        install_interaction(fh);
    end
    state = getappdata(fh, 'canlab_orthviews_state');

    need_load = isempty(state) || ~isfield(state, 'underlay') ...
                || isempty(state.underlay) ...
                || (~isempty(overlay_file) && ...
                    ~strcmp(state.underlay.file, resolve_underlay(overlay_file)));

    if need_load
        if isempty(overlay_file)
            file = canlab_get_underlay_image;
        else
            file = resolve_underlay(overlay_file);
        end
        u = load_underlay(file);
        if isempty(state)
            state = default_state(fh);
        end
        state.underlay = u;
        state.centre  = [0 0 0];
        state.zoom    = Inf;
        if ~isfield(state, 'blobs') || isempty(state.blobs)
            state.blobs = {};
        end
        if ~isfield(state, 'show_xhairs') || isempty(state.show_xhairs)
            state.show_xhairs = true;
        end
        setappdata(fh, 'canlab_orthviews_state', state);
        redraw_all(fh);
    end
end


function fh = find_canlab_orthviews_figure()
    fh = findobj('Type','figure','Tag','canlab_orthviews');
    if isempty(fh)
        fh = [];
    else
        fh = fh(1);
    end
end


function state = default_state(fh)
    state = struct( ...
        'fh', fh, ...
        'ax', [], ...
        'img', [], ...
        'xh_h', [], ...
        'xh_v', [], ...
        'centre', [0 0 0], ...
        'zoom', Inf, ...
        'show_xhairs', true, ...
        'show_colorbar', true, ...     % horizontal colorbar for autocolor blobs
        'bgcolor', [1 1 1], ...        % white background by default
        'underlay', [], ...
        'blobs', {{}}, ...
        'drag_panel', [], ...
        'drag_active', false );
end


function fh = build_figure()
    % 1x3 panels via create_figure; tag the figure for both our own lookup
    % and for legacy compatibility with code that searches for SPM's
    % 'Graphics' tag. Default background is white; can be overridden via
    % the 'backgroundcolor' / 'black' options.
    [fh, ax] = create_figure('canlab_orthviews', 1, 3);
    bg = [1 1 1];
    set(fh, 'Tag', 'canlab_orthviews', ...
            'Name', 'CANlab orthviews', ...
            'UserData', 'Graphics', ...
            'Color', bg);

    titles = {'Sagittal','Coronal','Axial'};
    w = 1/3;
    for k = 1:3
        axis(ax(k), 'image','off');
        set(ax(k), ...
            'Color', bg, ...
            'XColor','none','YColor','none', ...
            'Units','normalized', ...
            'LooseInset',[0 0 0 0]);
        title(ax(k), titles{k}, 'Color', fg_for(bg), 'FontSize', 11);
        hold(ax(k), 'on');

        % Tight per-axis position: no horizontal padding so panels are
        % flush; reserve a strip at the bottom of the figure for the
        % autocolor colorbar.
        set(ax(k), 'Position', [(k-1)*w, 0.11, w, 0.85]);

        % Hide the modern axes toolbar (and zoom/pan default interactions)
        % so clicks always reach our ButtonDownFcn cleanly.
        try, ax(k).Toolbar.Visible = 'off'; end %#ok<TRYNC>
        try, disableDefaultInteractivity(ax(k)); end %#ok<TRYNC>
    end

    % Initial state
    state = default_state(fh);
    state.ax = ax;
    state.bgcolor = bg;
    setappdata(fh, 'canlab_orthviews_state', state);

    % Click + right-click context menu wiring (attach to each axis)
    install_interaction(fh);
end


function rebuild_panels(fh)
    % Rebuild the three sagittal/coronal/axial panels inside an existing
    % canlab_orthviews figure (used to recover from external clobbering).
    state = getappdata(fh, 'canlab_orthviews_state');
    if isempty(state), return; end
    bg = [1 1 1];
    if isfield(state, 'bgcolor') && ~isempty(state.bgcolor)
        bg = state.bgcolor;
    end

    % Remove any leftover plain axes (the colorbar axis has its own tag
    % and is preserved).
    leftovers = findall(fh, 'Type', 'axes');
    for k = 1:numel(leftovers)
        if ~strcmp(get(leftovers(k), 'Tag'), 'canlab_orthviews_cbar')
            delete(leftovers(k));
        end
    end

    titles = {'Sagittal','Coronal','Axial'};
    w = 1/3;
    ax = gobjects(1, 3);
    for k = 1:3
        ax(k) = axes('Parent', fh, ...
            'Units','normalized', ...
            'Position', [(k-1)*w, 0.11, w, 0.85], ...
            'Color', bg, ...
            'XColor','none','YColor','none', ...
            'LooseInset',[0 0 0 0]);
        axis(ax(k), 'image','off');
        title(ax(k), titles{k}, 'Color', fg_for(bg), 'FontSize', 11);
        hold(ax(k), 'on');
        try, ax(k).Toolbar.Visible = 'off'; end %#ok<TRYNC>
        try, disableDefaultInteractivity(ax(k)); end %#ok<TRYNC>
    end
    state.ax = ax;
    setappdata(fh, 'canlab_orthviews_state', state);

    % Click + context menu must be re-wired to the fresh axes
    install_interaction(fh);
end


function apply_bgcolor(fh)
    % Recolor the figure, all axes, and titles to match state.bgcolor.
    state = getappdata(fh, 'canlab_orthviews_state');
    if isempty(state), return; end
    if ~isfield(state, 'bgcolor') || isempty(state.bgcolor)
        state.bgcolor = [1 1 1];
        setappdata(fh, 'canlab_orthviews_state', state);
    end
    bg = state.bgcolor;
    fg = fg_for(bg);
    set(fh, 'Color', bg);
    for k = 1:numel(state.ax)
        if ~isgraphics(state.ax(k)), continue; end
        set(state.ax(k), 'Color', bg);
        th = get(state.ax(k), 'Title');
        if isgraphics(th), set(th, 'Color', fg); end
    end
end


function install_interaction(fh)
    % Wire click + right-click context menu on each panel. Uses
    % non-anonymous handles to local functions and stashes per-axis state
    % (panel_id) via appdata, so callbacks survive across reruns / file
    % edits without going stale ("Unable to find function @(s,e)...").
    state = getappdata(fh, 'canlab_orthviews_state');
    if isempty(state) || ~isfield(state, 'ax') || isempty(state.ax)
        return
    end

    % Drop any previously attached context menu before replacing it
    old = findobj(fh, 'Type', 'uicontextmenu', 'Tag', 'canlab_orthviews_cmenu');
    if ~isempty(old), delete(old); end

    cmenu = uicontextmenu(fh, 'Tag', 'canlab_orthviews_cmenu');
    uimenu(cmenu, 'Label','Zoom in (20 mm)',  'Tag','canlab_zoom_20', ...
        'Callback', @ctx_menu_cb);
    uimenu(cmenu, 'Label','Zoom in (40 mm)',  'Tag','canlab_zoom_40', ...
        'Callback', @ctx_menu_cb);
    uimenu(cmenu, 'Label','Zoom in (80 mm)',  'Tag','canlab_zoom_80', ...
        'Callback', @ctx_menu_cb);
    uimenu(cmenu, 'Label','Reset zoom (full)','Tag','canlab_zoom_inf', ...
        'Callback', @ctx_menu_cb);
    uimenu(cmenu, 'Label','Toggle crosshairs','Tag','canlab_toggle_xhairs', ...
        'Separator','on','Callback', @ctx_menu_cb);
    uimenu(cmenu, 'Label','Re-center at origin','Tag','canlab_recenter', ...
        'Callback', @ctx_menu_cb);
    uimenu(cmenu, 'Label','Clear blobs','Tag','canlab_clear_blobs', ...
        'Callback', @ctx_menu_cb);

    for k = 1:numel(state.ax)
        if ~isgraphics(state.ax(k)), continue; end
        setappdata(state.ax(k), 'canlab_orthviews_panel_id', k);
        set(state.ax(k), ...
            'ButtonDownFcn', @axis_click_cb, ...
            'UIContextMenu', cmenu);
    end
end


function axis_click_cb(src, evt)
    % Single named entry point so MATLAB can re-resolve the local
    % function at click time. Reads the panel id from appdata, performs
    % the initial reposition, and then installs figure-level motion +
    % up callbacks so the user can drag the crosshair through slices.
    fh = ancestor(src, 'figure');
    if isempty(fh) || ~isgraphics(fh), return; end
    panel_id = getappdata(src, 'canlab_orthviews_panel_id');
    if isempty(panel_id), return; end

    % Initial click — render at full resolution
    on_axis_click(src, evt, fh, panel_id);

    % Begin drag: subsequent mouse-moves on the figure will reposition;
    % the renderer drops to npix=128 while drag_active is true.
    state = getappdata(fh, 'canlab_orthviews_state');
    state.drag_panel  = panel_id;
    state.drag_active = true;
    setappdata(fh, 'canlab_orthviews_state', state);

    set(fh, 'WindowButtonMotionFcn', @drag_motion_cb, ...
            'WindowButtonUpFcn',     @drag_up_cb);
end


function drag_motion_cb(src, ~)
    fh = ancestor(src, 'figure'); if isempty(fh), fh = src; end
    state = getappdata(fh, 'canlab_orthviews_state');
    if ~isfield(state, 'drag_panel') || isempty(state.drag_panel), return; end
    panel_id = state.drag_panel;
    if panel_id < 1 || panel_id > numel(state.ax), return; end
    ax = state.ax(panel_id);
    if ~isgraphics(ax), return; end

    cp = get(ax, 'CurrentPoint');
    u  = cp(1, 1);
    v  = cp(1, 2);

    xyz = state.centre;
    switch panel_id
        case 1, xyz(2) = u; xyz(3) = v;
        case 2, xyz(1) = u; xyz(3) = v;
        case 3, xyz(1) = u; xyz(2) = v;
    end
    state.centre = xyz;
    setappdata(fh, 'canlab_orthviews_state', state);
    redraw_all(fh);
end


function drag_up_cb(src, ~)
    fh = ancestor(src, 'figure'); if isempty(fh), fh = src; end
    set(fh, 'WindowButtonMotionFcn', '', 'WindowButtonUpFcn', '');
    state = getappdata(fh, 'canlab_orthviews_state');
    state.drag_panel  = [];
    state.drag_active = false;
    setappdata(fh, 'canlab_orthviews_state', state);
    % Re-render at full resolution now that the drag has ended
    redraw_all(fh);
end


function ctx_menu_cb(src, ~)
    % Single named callback for every context-menu item; dispatches on
    % the menu item's Tag.
    fh = ancestor(src, 'figure');
    if isempty(fh) || ~isgraphics(fh), return; end
    switch get(src, 'Tag')
        case 'canlab_zoom_20',       ctx_zoom(fh, 20);
        case 'canlab_zoom_40',       ctx_zoom(fh, 40);
        case 'canlab_zoom_80',       ctx_zoom(fh, 80);
        case 'canlab_zoom_inf',      ctx_zoom(fh, Inf);
        case 'canlab_toggle_xhairs', ctx_toggle_xhairs(fh);
        case 'canlab_recenter',      ctx_recenter(fh, [0 0 0]);
        case 'canlab_clear_blobs',   ctx_clear_blobs(fh);
    end
end


function ctx_zoom(fh, mm)
    state = getappdata(fh, 'canlab_orthviews_state');
    state.zoom = mm;
    setappdata(fh, 'canlab_orthviews_state', state);
    redraw_all(fh);
end

function ctx_toggle_xhairs(fh)
    state = getappdata(fh, 'canlab_orthviews_state');
    state.show_xhairs = ~state.show_xhairs;
    setappdata(fh, 'canlab_orthviews_state', state);
    redraw_crosshairs(fh);
end

function ctx_recenter(fh, xyz)
    state = getappdata(fh, 'canlab_orthviews_state');
    state.centre = xyz;
    setappdata(fh, 'canlab_orthviews_state', state);
    redraw_all(fh);
end

function ctx_clear_blobs(fh)
    state = getappdata(fh, 'canlab_orthviews_state');
    state.blobs = {};
    setappdata(fh, 'canlab_orthviews_state', state);
    redraw_all(fh);
end


function on_axis_click(src, evt, fh, panel_id)
    % src is the axes that was clicked. Compute the world (mm) coordinate
    % from the click position (or from CurrentPoint as a fallback) and
    % reposition the crosshair accordingly.
    state = getappdata(fh, 'canlab_orthviews_state');

    % Prefer event-data IntersectionPoint (modern HG2 hit events). This
    % also lets tests inject a synthetic point without writing to
    % CurrentPoint, which is read-only in newer MATLAB versions.
    [u, v] = click_point(src, evt);

    % Map (u, v) on the panel back to (x, y, z) mm using the current centre
    xyz = state.centre;
    switch panel_id
        case 1   % sagittal: panel axes are (Y mm, Z mm)
            xyz(2) = u;
            xyz(3) = v;
        case 2   % coronal: panel axes are (X mm, Z mm)
            xyz(1) = u;
            xyz(3) = v;
        case 3   % axial: panel axes are (X mm, Y mm)
            xyz(1) = u;
            xyz(2) = v;
    end
    state.centre = xyz;
    setappdata(fh, 'canlab_orthviews_state', state);
    redraw_all(fh);
end


function [u, v] = click_point(src, evt)
    % Resolve the (u, v) world coordinates of a click event on an axis.
    % Honors evt.IntersectionPoint when available (HG2), else falls back
    % to the axis CurrentPoint. Accepts a plain struct with the field as
    % a synthetic event for tests.
    u = []; v = [];
    if ~isempty(evt) && isstruct(evt) && isfield(evt, 'IntersectionPoint')
        p = evt.IntersectionPoint;
        u = p(1); v = p(2);
        return;
    end
    if ~isempty(evt) && isobject(evt) ...
            && isprop(evt, 'IntersectionPoint') ...
            && ~isempty(evt.IntersectionPoint)
        p = evt.IntersectionPoint;
        u = p(1); v = p(2);
        return;
    end
    cp = get(src, 'CurrentPoint');
    u = cp(1, 1);
    v = cp(1, 2);
end


function s = recenter_on_first_blob(s)
    % Recenter the crosshair on the most-prominent blob. With one combined
    % layer (autocolor / solid), this is just that layer's peak. With one
    % layer per cluster ('unique' / 'posneg'), pick the layer with the
    % largest support; for uniform-valued layers (typical of 'unique'),
    % use that cluster's centroid; otherwise its peak voxel.

    if isempty(s.blobs), return; end

    best = 1; best_n = -1;
    for k = 1:numel(s.blobs)
        L = s.blobs{k};
        if ~isfield(L, 'vol') || isempty(L.vol), continue; end
        n = nnz(L.vol ~= 0 & isfinite(L.vol));
        if n > best_n, best_n = n; best = k; end
    end

    L = s.blobs{best};
    if isfield(L, 'vol') && ~isempty(L.vol)
        Vabs = abs(L.vol);
        Vabs(~isfinite(Vabs)) = 0;
        nz = Vabs(Vabs > 0);
        uvals = unique(nz);
        if numel(uvals) <= 1
            % Uniform values across the cluster -> centroid is the only
            % meaningful focal point.
            mask = Vabs > 0;
            [I, J, K] = ind2sub(size(L.vol), find(mask));
            mm = L.M * [mean(I); mean(J); mean(K); 1];
        else
            [~, idx] = max(Vabs(:));
            [i, j, k] = ind2sub(size(L.vol), idx);
            mm = L.M * [i; j; k; 1];
        end
        s.centre = reshape(mm(1:3), 1, []);
        return
    end

    if isfield(L, 'XYZmm') && ~isempty(L.XYZmm)
        s.centre = reshape(mean(L.XYZmm, 2), 1, []);
    end
end


% =========================================================================
% Underlay I/O
% =========================================================================

function file = resolve_underlay(arg)
    if isempty(arg)
        file = canlab_get_underlay_image;
        return;
    end
    if ischar(arg) || isstring(arg)
        arg = char(arg);
        if exist(arg, 'file') == 2
            file = arg;
        else
            file = canlab_get_underlay_image(arg);
        end
    else
        file = canlab_get_underlay_image;
    end
end


function u = load_underlay(file)
    if isempty(file) || ~exist(file, 'file')
        error('canlab_orthviews:underlay', ...
            'Underlay image not found: %s', file);
    end

    [vol, M, dim] = read_nifti_volume(file);

    % Use a robust display window (2nd–98th percentile of non-zero voxels)
    finite_nz = vol(isfinite(vol) & vol ~= 0);
    if isempty(finite_nz)
        win = [0 1];
    else
        win = double(prctile(finite_nz, [2 99]));
        if win(1) == win(2), win = [min(finite_nz) max(finite_nz)]; end
        if win(1) == win(2), win = win + [0 1]; end
    end

    u = struct( ...
        'file',   file, ...
        'vol',    single(vol), ...     % 3-D underlay
        'mat',    M, ...               % 4x4 affine, 1-based, column-vec convention
        'dim',    dim, ...
        'window', win, ...
        'fovmm',  compute_fov(M, dim));
end


function fov = compute_fov(M, dim)
    % Bounding box in mm of all eight voxel-corner positions
    [I, J, K] = ndgrid([1 dim(1)], [1 dim(2)], [1 dim(3)]);
    P = [I(:)'; J(:)'; K(:)'; ones(1, 8)];
    Q = M * P;
    fov = [min(Q(1:3, :), [], 2) max(Q(1:3, :), [], 2)];  % [3x2] [lo hi]
end


function [vol, M, dim] = read_nifti_volume(file)
    % Read a NIfTI / NIfTI.gz volume using MATLAB built-ins; convert affine
    % from MATLAB's row-vector convention (Transform.T) to the SPM-style
    % 1-based column-vector M used throughout CANlab.
    try
        info = niftiinfo(file);
        vol  = niftiread(info);
    catch ME
        % Fall back: maybe an Analyze .img with a .hdr companion
        [p, n, e] = fileparts(file);
        alt = '';
        if strcmpi(e, '.img')
            alt = fullfile(p, [n '.nii']);
        elseif strcmpi(e, '.nii')
            alt = fullfile(p, [n '.img']);
        end
        if ~isempty(alt) && exist(alt, 'file')
            info = niftiinfo(alt);
            vol  = niftiread(info);
        else
            rethrow(ME);
        end
    end

    vol = double(vol);
    if ndims(vol) > 3
        vol = vol(:, :, :, 1);  % use first volume if 4-D
    end
    dim = size(vol);
    if numel(dim) < 3, dim(end+1:3) = 1; end

    % info.Transform.T is row-vector convention: [x y z 1] = [i j k 1]*T,
    % with 0-based voxel indices in the NIfTI standard. Convert to:
    %   [x y z 1]' = M_1based * [i j k 1]' for 1-based indices.
    T0 = info.Transform.T';      % column-vector convention, 0-based
    % Compose with the 0->1 shift so callers can use 1-based indices.
    shift0to1 = [eye(3), -ones(3, 1); 0 0 0 1];
    M = T0 * shift0to1;
end


% =========================================================================
% Slice math
% =========================================================================

function [img2d, X, Y] = sample_slice(u, plane, fixed_mm, lim_x, lim_y, npix)
    % Resample the underlay onto a regular mm grid in `plane` at the fixed
    % third-coordinate `fixed_mm`. Returns the 2-D image, plus the mm
    % vectors X and Y used as the panel axes.

    X = linspace(lim_x(1), lim_x(2), npix);
    Y = linspace(lim_y(1), lim_y(2), npix);
    [XX, YY] = meshgrid(X, Y);

    nP = numel(XX);
    switch plane
        case 'sagittal'   % YZ plane at fixed X
            mm = [repmat(fixed_mm, 1, nP); XX(:)'; YY(:)'];
        case 'coronal'    % XZ plane at fixed Y
            mm = [XX(:)'; repmat(fixed_mm, 1, nP); YY(:)'];
        case 'axial'      % XY plane at fixed Z
            mm = [XX(:)'; YY(:)'; repmat(fixed_mm, 1, nP)];
    end
    mm4 = [mm; ones(1, nP)];
    vox = u.mat \ mm4;
    vi = vox(1, :); vj = vox(2, :); vk = vox(3, :);

    % Trilinear interpolation in voxel space
    img1d = interpolate_volume(u.vol, vi, vj, vk);
    img2d = reshape(img1d, size(XX));
end


function vals = interpolate_volume(V, i, j, k)
    [nx, ny, nz] = size(V);
    in = i >= 1 & i <= nx & j >= 1 & j <= ny & k >= 1 & k <= nz;
    vals = nan(size(i));
    if ~any(in), return; end

    ii = i(in); jj = j(in); kk = k(in);
    i0 = floor(ii); j0 = floor(jj); k0 = floor(kk);
    i1 = min(i0 + 1, nx); j1 = min(j0 + 1, ny); k1 = min(k0 + 1, nz);
    di = ii - i0; dj = jj - j0; dk = kk - k0;

    % Clip lower bound voxels (they may equal nx etc, force >=1)
    i0 = max(i0, 1); j0 = max(j0, 1); k0 = max(k0, 1);

    s = @(a,b,c) V(sub2ind(size(V), a, b, c));
    v000 = s(i0,j0,k0); v100 = s(i1,j0,k0); v010 = s(i0,j1,k0);
    v110 = s(i1,j1,k0); v001 = s(i0,j0,k1); v101 = s(i1,j0,k1);
    v011 = s(i0,j1,k1); v111 = s(i1,j1,k1);

    w = ...
        v000 .* ((1-di).*(1-dj).*(1-dk)) + v100 .* (di.*(1-dj).*(1-dk)) + ...
        v010 .* ((1-di).*dj.*(1-dk))     + v110 .* (di.*dj.*(1-dk))     + ...
        v001 .* ((1-di).*(1-dj).*dk)     + v101 .* (di.*(1-dj).*dk)     + ...
        v011 .* ((1-di).*dj.*dk)         + v111 .* (di.*dj.*dk);

    vals(in) = w;
end


% =========================================================================
% Rendering
% =========================================================================

function redraw_all(fh)
    state = getappdata(fh, 'canlab_orthviews_state');
    if isempty(state) || isempty(state.underlay), return; end

    % Defensive: external code that takes over `gcf` (e.g., subplot or
    % clf inside plot/montage helpers) can delete our panels even when
    % the canlab_orthviews figure handle is reused. Detect that and
    % rebuild the three panels in place before drawing.
    if ~isfield(state, 'ax') || numel(state.ax) ~= 3 ...
            || ~all(isgraphics(state.ax))
        rebuild_panels(fh);
        state = getappdata(fh, 'canlab_orthviews_state');
    end

    apply_bgcolor(fh);
    state = getappdata(fh, 'canlab_orthviews_state');

    u = state.underlay;
    xyz = state.centre;

    % Per-panel axis limits (FOV in mm, clamped by zoom)
    if isfinite(state.zoom) && state.zoom > 0
        z = state.zoom;
        lim_sag = [xyz(2) - z, xyz(2) + z;  xyz(3) - z, xyz(3) + z];   % Y, Z
        lim_cor = [xyz(1) - z, xyz(1) + z;  xyz(3) - z, xyz(3) + z];   % X, Z
        lim_axi = [xyz(1) - z, xyz(1) + z;  xyz(2) - z, xyz(2) + z];   % X, Y
    else
        lim_sag = [u.fovmm(2,1) u.fovmm(2,2); u.fovmm(3,1) u.fovmm(3,2)];
        lim_cor = [u.fovmm(1,1) u.fovmm(1,2); u.fovmm(3,1) u.fovmm(3,2)];
        lim_axi = [u.fovmm(1,1) u.fovmm(1,2); u.fovmm(2,1) u.fovmm(2,2)];
    end

    % Drop sampling resolution during click-and-drag so mouse-tracking
    % stays responsive; restored automatically on mouse-up.
    if isfield(state, 'drag_active') && state.drag_active
        npix = 128;
    else
        npix = 256;
    end
    [I1, X1, Y1] = sample_slice(u, 'sagittal', xyz(1), lim_sag(1,:), lim_sag(2,:), npix);
    [I2, X2, Y2] = sample_slice(u, 'coronal',  xyz(2), lim_cor(1,:), lim_cor(2,:), npix);
    [I3, X3, Y3] = sample_slice(u, 'axial',    xyz(3), lim_axi(1,:), lim_axi(2,:), npix);

    plot_panel(state.ax(1), I1, X1, Y1, u.window, state.bgcolor, 1);
    plot_panel(state.ax(2), I2, X2, Y2, u.window, state.bgcolor, 2);
    plot_panel(state.ax(3), I3, X3, Y3, u.window, state.bgcolor, 3);

    % Use the underlay's binary footprint to clip blob alpha so that
    % pre-smoothed blob tails can't spill outside the brain coverage.
    underlay_alpha = { ...
        double(isfinite(I1) & I1 ~= 0), ...
        double(isfinite(I2) & I2 ~= 0), ...
        double(isfinite(I3) & I3 ~= 0)};

    % Overlay blobs
    overlay_blobs(state, lim_sag, lim_cor, lim_axi, npix, underlay_alpha);

    % Crosshairs + xyz readout
    draw_crosshairs(state, xyz, lim_sag, lim_cor, lim_axi);
    update_title(state.ax(1), 'Sagittal', xyz);
    update_title(state.ax(2), 'Coronal',  xyz);
    update_title(state.ax(3), 'Axial',    xyz);

    % Horizontal colorbar for autocolor blob layers (no-op when no
    % autocolor layer is present, e.g. for 'unique' / 'posneg' / 'color')
    draw_colorbar(fh, state);
end


function draw_colorbar(fh, state)
    % Compute the data range across all autocolor blob layers; if none
    % are present, remove any existing colorbar axis and return.
    vmin = Inf;
    vmax = -Inf;
    has_auto = false;
    for k = 1:numel(state.blobs)
        L = state.blobs{k};
        if ~isfield(L, 'style') || ~strcmp(L.style, 'autocolor'), continue; end
        if ~isfield(L, 'vol') || isempty(L.vol), continue; end
        V = L.vol;
        in = isfinite(V) & V ~= 0;
        if ~any(in(:)), continue; end
        vmin = min(vmin, double(min(V(in))));
        vmax = max(vmax, double(max(V(in))));
        has_auto = true;
    end

    cb_ax = findall(fh, 'Type','axes','Tag','canlab_orthviews_cbar');

    if ~has_auto || ~isfinite(vmin) || ~isfinite(vmax) || (vmin == 0 && vmax == 0)
        if ~isempty(cb_ax), delete(cb_ax); end
        return
    end
    if isfield(state, 'show_colorbar') && ~state.show_colorbar
        if ~isempty(cb_ax), delete(cb_ax); end
        return
    end

    % Degenerate case: a single value collapses XLim. Expand a hair so
    % MATLAB accepts the limits while still labeling the lone value.
    if vmax <= vmin
        if vmin == 0
            vmin = -eps; vmax = eps;
        else
            d = abs(vmax) * 0.05;
            if d == 0, d = eps; end
            vmax = vmin + d;
        end
    end

    % Reuse the colorbar axis if it already exists; otherwise create one.
    if isempty(cb_ax)
        cb_ax = axes('Parent', fh, 'Tag', 'canlab_orthviews_cbar', ...
                     'Units', 'normalized');
    else
        cb_ax = cb_ax(1);
        cla(cb_ax);
    end

    % Hot/cool gradient matching the blob renderer's mapping.
    n    = 256;
    zmax = max(abs(vmin), abs(vmax));
    if zmax <= 0, zmax = 1; end

    vals = linspace(vmin, vmax, n);
    R = zeros(1, n); G = R; B = R;
    pos = vals >= 0;
    neg = ~pos;
    if any(pos)
        t = min(vals(pos) / zmax, 1);
        R(pos) = 1;  G(pos) = max(0, 1 - t);  B(pos) = 0;
    end
    if any(neg)
        t = min(-vals(neg) / zmax, 1);
        R(neg) = 0;  G(neg) = max(0, 1 - t);  B(neg) = 1;
    end
    rgb = cat(3, R, G, B);

    set(cb_ax, 'Position', [0.35, 0.04, 0.30, 0.045], ...
               'Color', state.bgcolor, ...
               'XLim', [vmin vmax], 'YLim', [0 1]);

    image('Parent', cb_ax, ...
          'XData', [vmin vmax], 'YData', [0 1], 'CData', rgb, ...
          'HitTest', 'off', 'PickableParts', 'none');

    fg = fg_for(state.bgcolor);
    set(cb_ax, ...
        'XLim', [vmin vmax], 'YLim', [0 1], ...
        'YTick', [], ...
        'XTick', [vmin vmax], ...
        'XTickLabel', {sig1_str(vmin), sig1_str(vmax)}, ...
        'TickLength', [0 0], ...
        'XColor', fg, 'YColor', fg, ...
        'FontSize', 10, ...
        'Box', 'on', ...
        'TickDir', 'out');

    % The colorbar isn't interactive — keep it out of the click path.
    try, cb_ax.Toolbar.Visible = 'off'; end %#ok<TRYNC>
    try, disableDefaultInteractivity(cb_ax); end %#ok<TRYNC>
    set(cb_ax, 'HitTest', 'off', 'PickableParts', 'none');
end


function s = sig1_str(x)
    % Format a number to 1 significant digit. Handles zero, small
    % magnitudes (no scientific notation up to 1e4), and large
    % magnitudes (scientific notation beyond).
    if x == 0, s = '0'; return; end
    p = floor(log10(abs(x)));
    r = round(x / 10^p) * 10^p;
    if abs(r) >= 1 && abs(r) < 1e4
        % integer-ish: drop trailing decimals
        s = sprintf('%g', r);
    else
        s = sprintf('%.1g', r);
    end
end


function plot_panel(ax, I, X, Y, win, bg, panel_id)
    % Make sure click + context-menu wiring survives a redraw by drawing
    % into a child image but keeping the axes object stable. Zero / NaN
    % underlay pixels are made transparent so outside-brain regions
    % inherit the figure background color.
    cla(ax);

    A = double(isfinite(I) & I ~= 0);

    img = image('Parent', ax, ...
                'XData', X, 'YData', Y, 'CData', I, ...
                'AlphaData', A, ...
                'CDataMapping', 'scaled');
    set(ax, 'YDir', 'normal', 'CLim', win, ...
            'XLim', [X(1) X(end)], 'YLim', [Y(1) Y(end)], ...
            'Color', bg, ...
            'DataAspectRatio', [1 1 1]);
    set(ax, 'XTick', [], 'YTick', [], 'Box', 'off', ...
            'XColor', 'none', 'YColor', 'none');
    colormap(ax, gray(256));

    % Image should not steal clicks from the axis itself
    set(img, 'HitTest', 'off', 'PickableParts', 'none');
    set(ax,  'HitTest', 'on',  'PickableParts', 'all');

    title(ax, '', 'Color', fg_for(bg));
    setappdata(ax, 'canlab_orthviews_panel_id', panel_id);

    % Keep the axis toolbar / default interactions out of the way on
    % every redraw (cla can re-enable them).
    try, ax.Toolbar.Visible = 'off'; end %#ok<TRYNC>
    try, disableDefaultInteractivity(ax); end %#ok<TRYNC>
end


function overlay_blobs(state, lim_sag, lim_cor, lim_axi, npix, underlay_alpha)
    if isempty(state.blobs), return; end
    xyz = state.centre;
    if nargin < 6, underlay_alpha = {[], [], []}; end

    for b = 1:numel(state.blobs)
        L = state.blobs{b};
        if ~isfield(L, 'vol') || isempty(L.vol), continue; end

        plot_blob_panel(state.ax(1), 'sagittal', xyz(1), L, lim_sag, npix, underlay_alpha{1});
        plot_blob_panel(state.ax(2), 'coronal',  xyz(2), L, lim_cor, npix, underlay_alpha{2});
        plot_blob_panel(state.ax(3), 'axial',    xyz(3), L, lim_axi, npix, underlay_alpha{3});
    end
end


function plot_blob_panel(ax, plane, fixed_mm, L, limpair, npix, underlay_alpha)
    % Sample the blob volume on the same mm grid used for the underlay,
    % returning both interpolated values and a fractional 0..1 mask. When
    % L.smooth_edges is true, sampling is trilinear and the fractional
    % mask is used as alpha — this gives anti-aliased blob edges. When
    % false, sampling is nearest-neighbor and the alpha is binary, giving
    % crisp voxel boundaries.
    %
    % If `underlay_alpha` is provided (same size as the panel grid), the
    % blob alpha is multiplied by it, so pre-smoothed blob tails can't
    % render in the outside-brain region.

    if nargin < 7, underlay_alpha = []; end

    % Back-compat: layers built by older versions of this file may not
    % carry the smoothing fields. Default to smoothing on.
    if ~isfield(L, 'smooth_edges') || isempty(L.smooth_edges)
        L.smooth_edges = true;
    end

    [val2d, mask2d, X, Y] = sample_blob_slice(L, plane, fixed_mm, limpair, npix);
    if ~any(mask2d(:) > 0), return; end

    [rgb, A] = blob_slice_to_rgba(val2d, mask2d, L);

    if ~isempty(underlay_alpha) && isequal(size(underlay_alpha), size(A))
        A = A .* underlay_alpha;
    end

    image('Parent', ax, ...
          'XData', X, 'YData', Y, ...
          'CData', rgb, 'AlphaData', A, ...
          'HitTest', 'off', 'PickableParts', 'none');
end


function [val2d, mask2d, X, Y] = sample_blob_slice(L, plane, fixed_mm, limpair, npix)
    X = linspace(limpair(1,1), limpair(1,2), npix);
    Y = linspace(limpair(2,1), limpair(2,2), npix);
    [XX, YY] = meshgrid(X, Y);
    nP = numel(XX);

    switch plane
        case 'sagittal'   % YZ plane at fixed X
            mm = [repmat(fixed_mm, 1, nP); XX(:)'; YY(:)'];
        case 'coronal'    % XZ plane at fixed Y
            mm = [XX(:)'; repmat(fixed_mm, 1, nP); YY(:)'];
        case 'axial'      % XY plane at fixed Z
            mm = [XX(:)'; YY(:)'; repmat(fixed_mm, 1, nP)];
    end
    mm4 = [mm; ones(1, nP)];
    vox = L.M \ mm4;
    vi = vox(1, :); vj = vox(2, :); vk = vox(3, :);

    use_smooth = L.smooth_edges && ~strcmp(L.style, 'unique');

    % Layers built by older code may not have a separately-stored mask;
    % derive a binary one in that case.
    if isfield(L, 'mask') && ~isempty(L.mask)
        mask_vol = L.mask;
    else
        mask_vol = single(L.vol ~= 0 & isfinite(L.vol));
    end

    if use_smooth
        % Trilinear sample of the values (for color) and of the mask (for
        % alpha) — the mask was smoothed independently in make_blob_layer
        % when pre-blur is active, so the alpha falloff is a real
        % Gaussian profile and not the dilated binary tail of smooth3(V).
        val1d = interpolate_volume(L.vol, vi, vj, vk);
        val1d(~isfinite(val1d)) = 0;

        mask1d = interpolate_volume(mask_vol, vi, vj, vk);
        mask1d(~isfinite(mask1d)) = 0;
        mask1d = max(0, min(1, mask1d));
    else
        val1d  = sample_volume_nn(L.vol,     vi, vj, vk);
        mask1d = sample_volume_nn(mask_vol,  vi, vj, vk);
        mask1d = single(mask1d > 0);   % crisp boundary
    end

    val2d  = reshape(val1d,  size(XX));
    mask2d = reshape(mask1d, size(XX));
end


function vals = sample_volume_nn(V, i, j, k)
    [nx, ny, nz] = size(V);
    i = round(i); j = round(j); k = round(k);
    in = i >= 1 & i <= nx & j >= 1 & j <= ny & k >= 1 & k <= nz;
    vals = zeros(size(i), 'like', V);
    if any(in)
        vals(in) = V(sub2ind([nx ny nz], i(in), j(in), k(in)));
    end
end


function [rgb, A] = blob_slice_to_rgba(val2d, mask2d, L)
    [h, w] = size(val2d);
    R = zeros(h, w); G = R; B = R;
    A = zeros(h, w);

    in = mask2d > 0;

    switch L.style
        case {'solid','solidcolor'}
            c = L.color; if isempty(c), c = [1 1 0]; end
            R(in) = c(1); G(in) = c(2); B(in) = c(3);
            A(in) = double(mask2d(in)) * L.alpha;

        case 'unique'
            % NN-only path: each pixel carries an integer label. Build one
            % color per unique label.
            u = unique(val2d(in));
            try
                colors = scn_standard_colors(numel(u));
            catch
                colors = arrayfun(@(k) [rand rand rand], 1:numel(u), ...
                    'UniformOutput', false);
            end
            for kk = 1:numel(u)
                m = in & (val2d == u(kk));
                c = colors{kk};
                R(m) = c(1); G(m) = c(2); B(m) = c(3);
                A(m) = L.alpha;
            end

        otherwise   % 'autocolor' — hot/cool split centered on 0
            zmax = 1;
            if ~isempty(L.clim), zmax = L.clim(2); end
            if zmax <= 0, zmax = max(abs(val2d(in))); end
            if zmax <= 0, zmax = 1; end

            pos = in & val2d > 0;
            neg = in & val2d < 0;

            if any(pos(:))
                tpos = min(double(val2d(pos)) / zmax, 1);
                R(pos) = 1;
                G(pos) = max(0, 1 - tpos);
                B(pos) = 0;
            end
            if any(neg(:))
                tneg = min(double(-val2d(neg)) / zmax, 1);
                R(neg) = 0;
                G(neg) = max(0, 1 - tneg);
                B(neg) = 1;
            end
            A(in) = double(mask2d(in)) * L.alpha;
    end

    rgb = cat(3, R, G, B);
end


function draw_crosshairs(state, xyz, lim_sag, lim_cor, lim_axi)
    if ~state.show_xhairs, return; end
    col = [0 0.85 0];   % bright green, readable on light & dark backgrounds
    gap = 5;            % mm gap centered on the crosshair

    % Sagittal: panel axes are (Y, Z)
    add_split_h(state.ax(1), lim_sag(1,:), xyz(2), xyz(3), gap, col);
    add_split_v(state.ax(1), xyz(2), lim_sag(2,:), xyz(3), gap, col);
    % Coronal: panel axes are (X, Z)
    add_split_h(state.ax(2), lim_cor(1,:), xyz(1), xyz(3), gap, col);
    add_split_v(state.ax(2), xyz(1), lim_cor(2,:), xyz(3), gap, col);
    % Axial: panel axes are (X, Y)
    add_split_h(state.ax(3), lim_axi(1,:), xyz(1), xyz(2), gap, col);
    add_split_v(state.ax(3), xyz(1), lim_axi(2,:), xyz(2), gap, col);
end


function add_split_h(ax, xrange, cx, cy, gap, col)
    % Horizontal line at y = cy with a gap around x = cx
    if xrange(1) < cx - gap
        add_line(ax, [xrange(1), cx - gap], [cy, cy], col);
    end
    if cx + gap < xrange(2)
        add_line(ax, [cx + gap, xrange(2)], [cy, cy], col);
    end
end


function add_split_v(ax, cx, yrange, cy, gap, col)
    % Vertical line at x = cx with a gap around y = cy
    if yrange(1) < cy - gap
        add_line(ax, [cx, cx], [yrange(1), cy - gap], col);
    end
    if cy + gap < yrange(2)
        add_line(ax, [cx, cx], [cy + gap, yrange(2)], col);
    end
end


function add_line(ax, x, y, col)
    line('Parent', ax, ...
         'XData', x, 'YData', y, ...
         'Color', col, ...
         'LineWidth', 1, ...
         'HitTest', 'off', ...
         'PickableParts', 'none');
end


function redraw_crosshairs(fh)
    redraw_all(fh);
end


function update_title(ax, name, xyz)
    bg = get(ax, 'Color');
    title(ax, sprintf('%s  (%.0f, %.0f, %.0f)', name, xyz(1), xyz(2), xyz(3)), ...
        'Color', fg_for(bg), 'FontSize', 10);
end
