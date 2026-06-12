function htmlpath = canlab_niivue(obj, varargin)
% canlab_niivue Render a CANlab image object as an interactive NiiVue web page.
%
% Bridges a CanlabCore neuroimaging object (or a NIfTI file on disk) to a
% self-contained, browser-based NiiVue viewer. The object is treated as the
% OVERLAY (the data to colorize); an anatomical UNDERLAY is resolved
% automatically (or supplied). The function writes an .html page that boots
% the vendored NiiVue 0.57.0 library and the CANlab viewer module, then
% optionally opens it in your web browser.
%
% Two output modes are supported:
%   - 'standalone' (default): a SINGLE self-contained .html file with the
%     image data base64-embedded and the NiiVue library, viewer JS, and CSS
%     inlined. Portable; can be emailed or dropped on a web server as-is.
%   - folder bundle ('standalone', false): an .html file plus copies of
%     niivue.js, canlab_niivue_viewer.js, and canlab_niivue.css, with the
%     image data written to data/*.nii.gz and referenced by relative URL.
%
% :Usage:
% ::
%
%     htmlpath = canlab_niivue(obj, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2026  Tor Wager and CANlab
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
%   **obj:**
%        The OVERLAY image (the data to colorize). One of:
%          - an fmri_data, statistic_image, atlas, or other image_vector
%            object, or
%          - a char/string path to an existing .nii or .nii.gz file.
%
% :Optional Inputs:
%
%   **'underlay', [char path or object]:**
%        Anatomical underlay. Default [] -> resolved via
%        canlab_get_underlay_image (the lab default template). May be a
%        filename, a known underlay keyword, or an image object (which is
%        written to a temp .nii.gz). Pass '' or 'none' to skip the underlay
%        and show the overlay alone.
%
%   **'outdir', [char]:**
%        Output directory for the html page and (folder mode) assets.
%        Default fullfile(pwd, 'canlab_niivue_output'). Created if missing.
%
%   **'fname', [char]:**
%        Name of the html file to write. Default 'index.html'.
%
%   **'title', [char]:**
%        Page title and visible H1 header. Default 'CANlab NiiVue viewer'.
%
%   **'standalone', [logical]:**
%        true (default) -> single self-contained html with base64-embedded
%        data and inlined niivue.js + viewer + css. false -> folder bundle
%        with relative-URL imports and data/*.nii.gz files.
%
%   **'colormap', [char]:**
%        NiiVue colormap for the OVERLAY positive values. Default 'inferno'.
%
%   **'colormapNegative', [char]:**
%        NiiVue colormap for the OVERLAY negative values (the CANlab
%        hot/cool split). Default 'winter'.
%
%   **'cal_min', [scalar]:**
%        Lower display threshold for the overlay AND the transparency cutoff:
%        voxels with |value| < cal_min render transparent so the underlay shows
%        through. Default [] -> auto from the overlay's smallest suprathreshold
%        magnitude (for a thresholded statistic_image, that is the threshold).
%
%   **'cal_max', [scalar]:**
%        Upper display limit for the overlay colormap. Default [] -> auto from
%        the overlay's largest suprathreshold magnitude.
%
%   **'opacity', [scalar]:**
%        Overlay opacity in [0 1]. Default 0.8 (semi-transparent, so anatomy
%        reads through the blobs). The page also has a live opacity slider.
%
%   **'underlay_2mm', [logical]:**
%        In standalone mode, resample a finer-than-2mm underlay down to ~2mm so
%        the single self-contained .html stays small (a 1mm template would add
%        ~4MB of base64). Default true. Ignored in folder mode (the underlay is
%        a separate cached file there, so full resolution is cheap).
%
%   **'thresh', [logical]:**
%        For statistic_image overlays, write the THRESHOLDED map (zeroing
%        non-significant voxels) rather than the raw values. Default true
%        for statistic_image, false otherwise. Ignored for non-objects.
%
%   **'open', [logical]:**
%        Open the resulting page in the web browser (web(htmlpath)).
%        Default true. 'noopen' to turn off.
%
%   **'doplot', [logical flag]:**
%        Synonym for 'open'; render/open the page. Default = true.
%        'noplot' to turn off.
%
%   **'doverbose', [logical flag]:**
%        Verbose output; default = true. 'noverbose' to turn off.
%
% :Outputs:
%
%   **htmlpath:**
%        Full path to the written .html page.
%
% :Examples:
% ::
%
%    % Threshold a t-map and view it in the browser:
%    dat = load_image_set('emotionreg');
%    t = threshold(ttest(dat), .005, 'unc', 'k', 10);
%    canlab_niivue(t)
%
%    % A folder bundle with a custom underlay and fixed color limits:
%    canlab_niivue(t, 'standalone', false, 'underlay', 'mni152_1mm', ...
%        'cal_min', 2, 'cal_max', 6, 'outdir', '~/niivue_demo');
%
%    % Render a NIfTI file directly, without opening a browser:
%    htmlpath = canlab_niivue('my_stat_map.nii.gz', 'noopen');
%
% :References:
%   NiiVue: Hanayik et al., niivue.com (BSD-2-Clause). CANlab tools:
%   https://canlab.github.io
%
% :See also:
%   - canlab_get_underlay_image
%   - canlab_results_fmridisplay
%   - fmridisplay, image_vector.surface, image_vector.montage

% ..
%    Programmers' notes:
%    2026/06  Created (Tor Wager / CANlab). Standalone mode base64-embeds
%             the .nii.gz bytes and inlines the vendored niivue.js plus the
%             viewer module (with its leading `import` line stripped, since
%             Niivue is already in scope when inlined).
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

% Parse special command keywords and remove them before inputParser

doplot = true;
plot_idx = strcmpi(varargin, 'plot');
if any(plot_idx)
    doplot = true;
    varargin(plot_idx) = [];
end
noplot_idx = strcmpi(varargin, 'noplot');
if any(noplot_idx)
    doplot = false;
    varargin(noplot_idx) = [];
end

% 'open' is a synonym for doplot here (open the page in a browser)
noopen_idx = strcmpi(varargin, 'noopen');
if any(noopen_idx)
    doplot = false;
    varargin(noopen_idx) = [];
end

verbose = true;
verbose_idx = strcmpi(varargin, 'noverbose');
if any(verbose_idx)
    verbose = false;
    varargin(verbose_idx) = [];   % remove so inputParser doesn't see it
end

% Default for 'thresh' depends on object type: statistic_image -> true.
default_thresh = isa(obj, 'statistic_image');

% Use inputParser to parse key/value pairs
% First add obligatory/non-conditional keywords

p = inputParser;
p.addRequired('obj', @(x) isa(x, 'image_vector') || ischar(x) || isstring(x));
p.addParameter('underlay', [], @(x) isempty(x) || ischar(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('outdir', fullfile(pwd, 'canlab_niivue_output'), @(x) ischar(x) || isstring(x));
p.addParameter('fname', 'index.html', @(x) ischar(x) || isstring(x));
p.addParameter('title', 'CANlab NiiVue viewer', @(x) ischar(x) || isstring(x));
p.addParameter('standalone', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('colormap', 'inferno', @(x) ischar(x) || isstring(x));
p.addParameter('colormapNegative', 'winter', @(x) ischar(x) || isstring(x));
p.addParameter('cal_min', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('cal_max', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('opacity', 0.8, @(x) isnumeric(x) && isscalar(x));
p.addParameter('underlay_2mm', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('thresh', default_thresh, @(x) islogical(x) || isnumeric(x));
p.addParameter('open', doplot, @(x) islogical(x) || isnumeric(x));

% Special key/value pairs that we have potentially set with optional keywords
p.addParameter('doplot', doplot, @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose', verbose, @(x) islogical(x) || isnumeric(x));

% process inputs and deal out to variables in workspace
p.parse(obj, varargin{:});

ARGS = p.Results;

% Get all field names in ARGS
fn = fieldnames(ARGS);

% Loop over fields and assign variables in THIS workspace
for i = 1:numel(fn)
    eval(sprintf('%s = ARGS.(''%s'');', fn{i}, fn{i}));
end

% END INPUT PARSER TEMPLATE

% Reconcile doplot / open: if either explicitly given, 'open' wins as the
% browser-launch control. doplot defaults to the same value. Read 'open' from
% ARGS (not the eval-created var) since it shadows the open() builtin.
doopen = logical(ARGS.open) && logical(doplot);

% Normalize string args to char. Read straight from ARGS (rather than the
% eval-created workspace vars) so the code analyzer can see the definitions
% and to avoid shadowing the title()/colormap() builtins.
outdir   = char(ARGS.outdir);
fname    = char(ARGS.fname);
ttl      = char(ARGS.title);
cmap     = char(ARGS.colormap);
cmapNeg  = char(ARGS.colormapNegative);

% Numeric overlay options (read from ARGS so the analyzer sees the definitions;
% cal_min/cal_max may be [] here and are auto-filled below).
cal_min  = ARGS.cal_min;
cal_max  = ARGS.cal_max;
opacity  = ARGS.opacity;
dolite   = logical(ARGS.underlay_2mm);

% -------------------------------------------------------------------------
% Locate asset directory
% -------------------------------------------------------------------------

thisdir  = fileparts(mfilename('fullpath'));
assetdir = fullfile(thisdir, 'assets');

template_file = fullfile(assetdir, 'canlab_niivue_template.html');
viewer_file   = fullfile(assetdir, 'canlab_niivue_viewer.js');
css_file      = fullfile(assetdir, 'canlab_niivue.css');
niivue_file   = fullfile(assetdir, 'niivue.js');

assert(exist(template_file, 'file') == 2, 'canlab_niivue: missing asset %s', template_file);
assert(exist(viewer_file, 'file') == 2, 'canlab_niivue: missing asset %s', viewer_file);
assert(exist(css_file, 'file') == 2, 'canlab_niivue: missing asset %s', css_file);
assert(exist(niivue_file, 'file') == 2, 'canlab_niivue: missing vendored library %s', niivue_file);

% -------------------------------------------------------------------------
% Prepare output directory
% -------------------------------------------------------------------------

if ~exist(outdir, 'dir'), mkdir(outdir); end
htmlpath = fullfile(outdir, fname);

% Temp dir for any objects we have to write out before encoding
tmpdir = tempname;
mkdir(tmpdir);
cleanupTmp = onCleanup(@() rmtmp(tmpdir));

% -------------------------------------------------------------------------
% Resolve OVERLAY to a .nii.gz file on disk
% -------------------------------------------------------------------------

dowrite_thresh = logical(thresh);

% Auto color limits / transparency cutoff from the overlay's suprathreshold
% magnitude, unless the caller specified them. cal_min doubles as the
% transparency threshold in the viewer (|value| < cal_min renders transparent).
if isa(obj, 'image_vector') && (isempty(cal_min) || isempty(cal_max))
    [cmin_auto, cmax_auto] = local_auto_callimits(obj, dowrite_thresh);
    if isempty(cal_min), cal_min = cmin_auto; end
    if isempty(cal_max), cal_max = cmax_auto; end
end

overlay_file = local_resolve_image(obj, fullfile(tmpdir, 'overlay.nii.gz'), ...
    dowrite_thresh, verbose);

% -------------------------------------------------------------------------
% Resolve UNDERLAY to a file path (or [] for none)
% -------------------------------------------------------------------------

has_underlay = true;

if isa(underlay, 'image_vector')
    underlay_file = local_resolve_image(underlay, fullfile(tmpdir, 'underlay.nii.gz'), false, verbose);

elseif (ischar(underlay) || isstring(underlay)) && ...
        (strcmpi(underlay, 'none') || strcmp(underlay, ''))
    has_underlay = false;
    underlay_file = '';

elseif ischar(underlay) || isstring(underlay)
    % Keyword or filename -> resolve via canlab_get_underlay_image
    underlay_file = canlab_get_underlay_image(char(underlay));
    if isempty(underlay_file) || ~exist(underlay_file, 'file')
        warning('canlab_niivue: could not resolve underlay; showing overlay alone.');
        has_underlay = false;
        underlay_file = '';
    end

else
    % Default [] -> lab default underlay
    underlay_file = canlab_get_underlay_image;
    if isempty(underlay_file) || ~exist(underlay_file, 'file')
        warning('canlab_niivue: could not resolve default underlay; showing overlay alone.');
        has_underlay = false;
        underlay_file = '';
    end
end

% For standalone, downsample a finer-than-2mm underlay to ~2mm so the single
% embedded .html stays small. Folder mode keeps the full-resolution file.
if has_underlay && logical(standalone) && dolite
    lite = local_lite_underlay(underlay_file, verbose);
    if isa(lite, 'image_vector')
        underlay_file = local_resolve_image(lite, fullfile(tmpdir, 'underlay.nii.gz'), false, verbose);
    end
end

% -------------------------------------------------------------------------
% Read template, viewer, css
% -------------------------------------------------------------------------

template_src = fileread(template_file);
viewer_src   = fileread(viewer_file);
css_src      = fileread(css_file);

% Strip the leading `import { ... } from './niivue.js';` line from the
% viewer for standalone inlining (Niivue is already in scope there). The
% contract guarantees this is the verbatim first line of the file.
viewer_inline = local_strip_first_import(viewer_src);

% -------------------------------------------------------------------------
% Build a small JS object literal for the overlay color/threshold options.
% Shared by both modes.
% -------------------------------------------------------------------------

opt_fields = sprintf('colormap: %s, colormapNegative: %s, opacity: %s', ...
    local_jsstr(cmap), local_jsstr(cmapNeg), num2str(opacity));

if ~isempty(cal_min)
    opt_fields = sprintf('%s, cal_min: %s', opt_fields, num2str(cal_min));
end
if ~isempty(cal_max)
    opt_fields = sprintf('%s, cal_max: %s', opt_fields, num2str(cal_max));
end

% -------------------------------------------------------------------------
% Assemble STYLE_BLOCK and SCRIPT_BLOCK per mode
% -------------------------------------------------------------------------

if logical(standalone)
    % ---------------- STANDALONE ----------------
    style_block = sprintf('<style>\n%s\n</style>', css_src);

    overlay_b64 = local_base64_of_file(overlay_file);

    if has_underlay
        underlay_b64 = local_base64_of_file(underlay_file);
        cfg = sprintf(['{\n', ...
            '  underlay: { base64: UNDERLAY_B64, name: "underlay.nii.gz" },\n', ...
            '  overlay:  { base64: OVERLAY_B64, name: "overlay.nii.gz" },\n', ...
            '  %s\n}'], opt_fields);
    else
        underlay_b64 = '';
        cfg = sprintf(['{\n', ...
            '  overlay: { base64: OVERLAY_B64, name: "overlay.nii.gz" },\n', ...
            '  %s\n}'], opt_fields);
    end

    % Embed base64 payloads as JS const string literals, then reference them.
    boot = [ ...
        sprintf('const OVERLAY_B64 = %s;\n', local_jsstr(overlay_b64))];
    if has_underlay
        boot = [boot sprintf('const UNDERLAY_B64 = %s;\n', local_jsstr(underlay_b64))];
    end
    boot = [boot sprintf('canlabNiivue("gl", %s);\n', cfg)];

    script_block = [ ...
        sprintf('<script type="module">\n') ...
        niivue_inline_marker() ...
        fileread(niivue_file) newline ...
        viewer_inline newline ...
        boot ...
        sprintf('</script>')];

else
    % ---------------- FOLDER BUNDLE ----------------
    style_block = '<link rel="stylesheet" href="canlab_niivue.css">';

    % Copy the three asset files into outdir
    copyfile(niivue_file, fullfile(outdir, 'niivue.js'));
    copyfile(viewer_file, fullfile(outdir, 'canlab_niivue_viewer.js'));
    copyfile(css_file,    fullfile(outdir, 'canlab_niivue.css'));

    % Write data files
    datadir = fullfile(outdir, 'data');
    if ~exist(datadir, 'dir'), mkdir(datadir); end
    copyfile(overlay_file, fullfile(datadir, 'overlay.nii.gz'));

    if has_underlay
        copyfile(underlay_file, fullfile(datadir, 'underlay.nii.gz'));
        cfg = sprintf(['{\n', ...
            '  underlay: "./data/underlay.nii.gz",\n', ...
            '  overlay:  "./data/overlay.nii.gz",\n', ...
            '  %s\n}'], opt_fields);
    else
        cfg = sprintf(['{\n', ...
            '  overlay: "./data/overlay.nii.gz",\n', ...
            '  %s\n}'], opt_fields);
    end

    script_block = sprintf([ ...
        '<script type="module">\n', ...
        'import { canlabNiivue } from "./canlab_niivue_viewer.js";\n', ...
        'canlabNiivue("gl", %s);\n', ...
        '</script>'], cfg);
end

% -------------------------------------------------------------------------
% Token substitution and write
% -------------------------------------------------------------------------

html = template_src;
html = strrep(html, '{{TITLE}}', ttl);
html = strrep(html, '{{STYLE_BLOCK}}', style_block);
html = strrep(html, '{{SCRIPT_BLOCK}}', script_block);

fid = fopen(htmlpath, 'w');
if fid == -1, error('canlab_niivue: could not open %s for writing.', htmlpath); end
fwrite(fid, html, 'char');
fclose(fid);

% -------------------------------------------------------------------------
% Report and open
% -------------------------------------------------------------------------

if verbose
    dashes = '----------------------------------------------';
    fprintf('%s\n', dashes);
    fprintf('canlab_niivue: wrote %s viewer\n', ternary(logical(standalone), 'standalone', 'folder-bundle'));
    fprintf('  HTML:     %s\n', htmlpath);
    if ~logical(standalone)
        fprintf('  Assets:   niivue.js, canlab_niivue_viewer.js, canlab_niivue.css\n');
        fprintf('  Data:     data/overlay.nii.gz%s\n', ternary(has_underlay, ', data/underlay.nii.gz', ''));
    end
    if has_underlay
        fprintf('  Underlay: %s\n', underlay_file);
    else
        fprintf('  Underlay: (none)\n');
    end
    fprintf('%s\n', dashes);
end

if doopen
    web(htmlpath, '-browser');
end

end % main function


% =========================================================================
% Subfunctions
% =========================================================================

function out_file = local_resolve_image(im, target_nii, do_thresh, verbose)
% Resolve an overlay/underlay argument to a gzipped .nii.gz file on disk.
% target_nii must end in '.nii.gz'; the returned file is always valid gzipped
% NIfTI, so its bytes match the '.nii.gz' name NiiVue uses for format detection.
%   - char/string path: if already '.nii.gz', returned as-is; if uncompressed
%     ('.nii'/'.img'), gzipped into target_nii.
%   - image_vector object: written to a plain '.nii' via write() (SPM's writer
%     rejects a '.gz' extension), then gzipped to target_nii.

if ischar(im) || isstring(im)
    in_file = char(im);
    if ~exist(in_file, 'file')
        error('canlab_niivue: image file not found: %s', in_file);
    end

    [~, ~, ext] = fileparts(in_file);
    if strcmpi(ext, '.gz')
        % Already gzipped (.nii.gz) -> use directly.
        out_file = in_file;
    else
        % Uncompressed NIfTI -> gzip into the target so the advertised
        % '.nii.gz' name matches the actual (gzipped) bytes.
        local_gzip_to_target(in_file, target_nii);
        out_file = target_nii;
    end
    return
end

% Object: write to a plain '.nii' (SPM rejects the '.gz' extension), then
% gzip to the '.nii.gz' target. statistic_image supports 'thresh'.
nii_tmp = strrep(target_nii, '.nii.gz', '.nii');

write_args = {'fname', nii_tmp, 'overwrite'};
if do_thresh && isa(im, 'statistic_image')
    write_args = [write_args, {'thresh'}];
end

if verbose
    fprintf('canlab_niivue: writing image to %s\n', target_nii);
end

write(im, write_args{:});

if ~exist(nii_tmp, 'file')
    error('canlab_niivue: write() did not produce expected file: %s', nii_tmp);
end

local_gzip_to_target(nii_tmp, target_nii);

out_file = target_nii;

if ~exist(out_file, 'file')
    error('canlab_niivue: failed to produce %s', out_file);
end

end % local_resolve_image


function local_gzip_to_target(src_nii, target_gz)
% Gzip src_nii so the result lands exactly at target_gz (which ends in .nii.gz).
% gzip() writes <basename>.gz into a folder, so we arrange the uncompressed
% name to be the target minus '.gz', gzip it, and clean up any temp copy.

[tdir, tname] = fileparts(target_gz);     % tname e.g. 'overlay.nii', drops '.gz'
plain_target  = fullfile(tdir, tname);    % e.g. <tdir>/overlay.nii

made_copy = ~strcmp(src_nii, plain_target);
if made_copy
    copyfile(src_nii, plain_target);
end

gzip(plain_target);                       % writes plain_target + '.gz' == target_gz

if made_copy && exist(plain_target, 'file')
    delete(plain_target);                 % remove the intermediate uncompressed copy
end

end % local_gzip_to_target


function out = local_lite_underlay(underlay_file, verbose)
% Resample a finer-than-2mm underlay down to ~2mm (keeps standalone HTML small).
% Returns an fmri_data object when resampled, or [] when left as-is (caller then
% keeps the original file). Uses a canonical 2mm template as the target grid.

out = [];

ref = '';
for c = {'avg152T1.nii', 'single_subj_T1.nii'}
    w = which(c{1});
    if ~isempty(w), ref = w; break; end
end
if isempty(ref), return; end   % no 2mm reference available -> keep original

try
    u = fmri_data(underlay_file);
    vs = abs([u.volInfo.mat(1, 1), u.volInfo.mat(2, 2), u.volInfo.mat(3, 3)]);
    if max(vs) >= 1.9
        return    % already ~2mm or coarser
    end
    if verbose
        fprintf('canlab_niivue: resampling underlay (%.2gmm) to ~2mm for standalone size.\n', max(vs));
    end
    out = resample_space(u, fmri_data(ref));
catch
    out = [];
end

end % local_lite_underlay


function [cmin, cmax] = local_auto_callimits(im, do_thresh)
% Robust overlay color limits from the suprathreshold magnitude of the data.
% cmin = smallest nonzero |value| (the effective threshold), cmax = largest.
% Returns [] if it cannot compute (caller then falls back to NiiVue auto).

cmin = [];
cmax = [];

try
    d = im.dat(:, 1);

    % For a thresholded statistic_image, zero out non-significant voxels so the
    % limits reflect only the blobs that will be displayed.
    if do_thresh && isa(im, 'statistic_image') && ~isempty(im.sig)
        s = logical(im.sig(:, 1));
        d(~s) = 0;
    end

    a = abs(d(:));
    a = a(a > 0 & ~isnan(a) & ~isinf(a));
    if isempty(a), return; end

    cmin = double(min(a));
    cmax = double(max(a));
    if cmax <= cmin
        cmax = cmin + max(eps(cmin), 1e-6);
    end
catch
    cmin = [];
    cmax = [];
end

end % local_auto_callimits


function s = local_base64_of_file(fpath)
% Base64-encode the raw bytes of a file.

fid = fopen(fpath, 'r');
if fid == -1, error('canlab_niivue: cannot read %s', fpath); end
bytes = fread(fid, Inf, '*uint8');
fclose(fid);

if exist('matlab.net.base64encode', 'file') || ...
        ~isempty(which('matlab.net.base64encode'))
    s = matlab.net.base64encode(bytes);
else
    % Fallback for older MATLAB: use Java
    s = char(org.apache.commons.codec.binary.Base64.encodeBase64String(bytes));
end

end % local_base64_of_file


function viewer_out = local_strip_first_import(viewer_src)
% Remove the first line of the viewer source (the niivue.js import).
% The contract guarantees the first line is exactly the import statement.

nl = find(viewer_src == newline, 1, 'first');
if isempty(nl)
    viewer_out = viewer_src;
else
    viewer_out = viewer_src(nl+1:end);
end

end % local_strip_first_import


function s = local_jsstr(str)
% Encode a MATLAB char array as a safe single-quoted JS string literal.

str = char(str);
str = strrep(str, '\', '\\');
str = strrep(str, '''', '\''');
str = strrep(str, newline, '\n');
str = strrep(str, sprintf('\r'), '\r');
s = ['''' str ''''];

end % local_jsstr


function m = niivue_inline_marker()
% A short comment marking the start of the inlined NiiVue bundle.
m = sprintf('/* --- inlined niivue.js (vendored 0.57.0) --- */\n');
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


function rmtmp(d)
try
    if exist(d, 'dir'), rmdir(d, 's'); end
catch
end
end
