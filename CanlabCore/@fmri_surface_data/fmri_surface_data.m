classdef fmri_surface_data < image_vector
% fmri_surface_data: CANlab data object for cortical-surface and grayordinate (CIFTI) data
%
% fmri_surface_data is the surface/grayordinate analogue of fmri_data. It wraps
% the HCP CIFTI grayordinate standard (cortical-surface vertices + subcortical
% voxels) as a true CANlab object that is interoperable with the rest of the
% toolbox, while holding data in the canonical flat [grayordinates x maps] .dat
% matrix so generic analysis methods (ica, descriptives, ...) work unchanged.
%
% It is a subclass of image_vector (NOT fmri_data): it inherits the .dat-based
% machinery but does NOT reproduce fmri_data's volInfo-centric, empty-voxel
% squeezing design. CIFTI data is already compact (the medial wall is excluded
% by construction and there are no out-of-brain voxels), so:
%   - .dat is ALWAYS the full [nGrayordinates x nMaps] matrix, in 1:1 row
%     correspondence with .brain_model. Rows are never dropped.
%   - remove_empty / replace_empty are no-ops (see those methods).
%   - removed_voxels / removed_images are vestigial all-false vectors.
% See docs/fmri_surface_data_design_plan.md (decision D5b).
%
% -------------------------------------------------------------------------
% Construction
% -------------------------------------------------------------------------
%   obj = fmri_surface_data                     % empty object
%   obj = fmri_surface_data(cifti_filename)     % .dscalar/.dtseries/.dlabel.nii
%   obj = fmri_surface_data(gifti_filename)     % .func/.shape/.label/.surf.gii
%   obj = fmri_surface_data(cifti_struct)       % output of canlab_read_cifti
%   obj = fmri_surface_data('dat', X, 'brain_model', bm, ...)   % key/value pairs
%   obj = fmri_surface_data(image_vector_obj)   % recast a matching object
%
% Native I/O is used (canlab_read_cifti / canlab_read_gifti) -- no external
% toolbox (gifti, FieldTrip, cifti-matlab, Connectome Workbench) is required.
%
% -------------------------------------------------------------------------
% Key properties
% -------------------------------------------------------------------------
%   dat            [nGrayordinates x nMaps] single  (inherited) data matrix.
%   brain_model    struct: the geometry source of truth (mirrors CIFTI
%                  BrainModels). Fields .models{i} (.struct, .type 'surf'|'vox',
%                  .start, .count, .numvert, .vertlist 0-based, .voxlist 3xN),
%                  .vol (.dims, .sform), .grayordinate_type, .cluster.
%   geom           struct: cortical mesh(es) for rendering/area (lazy).
%   imagetype      'dscalar'|'dtseries'|'dlabel'|'func'|'shape'|'label'.
%   series_info    struct for .dtseries (start/step/unit/exponent).
%   label_table    MATLAB table (key, name, rgba) for .dlabel/.label.
%   surface_space  e.g. 'fsLR_32k', 'fsaverage_164k'.
%   X / Y / covariates / images_per_session / metadata_table ...
%                  per-map annotations (same names/roles as fmri_data).
%
% :See also: fmri_data, image_vector, canlab_read_cifti, canlab_read_gifti

% ..
%    Author: CANlab, 2026. Part of the fmri_surface_data surface-data toolset.
%    See docs/fmri_surface_data_design_plan.md.
% ..

    properties

        % Geometry source of truth (mirrors the CIFTI BrainModels structure).
        % Replaces volInfo's role for surface + grayordinate data. See D3/D5b.
        brain_model = [];

        % Cortical mesh cache (faces/vertices per hemisphere) for rendering and
        % surface-area computations. Loaded lazily from bundled assets.
        geom = [];

        % CIFTI/GIFTI image type (a.k.a. CIFTI intent): dscalar | dtseries | dlabel | func | shape | label
        imagetype = '';

        % For .dtseries: struct with .start/.step/.unit/.exponent
        series_info = [];

        % For .dlabel/.label: MATLAB table with variables key, name, rgba
        label_table = [];

        % Canonical surface-space tag, e.g. 'fsLR_32k', 'fsaverage_164k'
        surface_space = '';

        % Optional same-length logical over grayordinates (or another
        % fmri_surface_data on the same space). Lightweight: masking just
        % zeros/selects rows -- no fmri_mask_image, no resampling (D5b).
        mask = [];

        % --- Per-map (column) annotations: same names/roles as fmri_data ---
        X                                   % design / predictor matrix [nMaps x p]
        Y = [];                             % outcomes [nMaps x q]
        Y_names;
        covariates;
        covariate_names = {''};
        images_per_session;
        metadata_table;
        image_metadata;
        additional_info = struct();         % free-form; also stashes source CIFTI xml/hdr

    end % properties

    methods

        function obj = fmri_surface_data(varargin)

            % Initialize defaults via the image_vector constructor
            obj = obj@image_vector();
            obj.removed_voxels = false(0, 1);
            obj.removed_images = false(0, 1);

            if nargin == 0, return; end

            % ---- Single positional argument: file / struct / object ----
            if nargin == 1
                a = varargin{1};

                if ischar(a) || isstring(a)
                    obj = load_surface_file(obj, char(a));
                    return

                elseif isstruct(a)
                    if isfield(a, 'cdata') && isfield(a, 'diminfo')
                        obj = from_cifti_struct(obj, a);     % canlab_read_cifti output
                    elseif isfield(a, 'vertices') || isfield(a, 'faces') || isfield(a, 'cdata')
                        obj = from_gifti_struct(obj, a);     % canlab_read_gifti output
                    else
                        obj = copy_matching_fields(obj, a);  % plain struct of properties
                    end
                    return

                elseif isa(a, 'image_vector')
                    obj = recast_from_object(obj, a);
                    return
                end
            end

            % ---- key / value pairs ----
            obj = parse_keyvalue(obj, varargin);

            % If a filename was passed via 'fullpath' and we have no data, load it
            if isempty(obj.dat) && ~isempty(obj.fullpath) && ischar(obj.fullpath) ...
                    && exist(obj.fullpath, 'file')
                obj = load_surface_file(obj, obj.fullpath);
            end

        end % constructor

    end % methods

    % ---------------------------------------------------------------------
    % Volume-only image_vector methods that have no meaningful surface /
    % grayordinate behavior. They are overridden here inside a (Hidden) block
    % so that (a) they no longer clutter methods(obj) / tab-completion, and
    % (b) calling one gives a clear message (pointing to the surface
    % equivalent) instead of a cryptic volume-machinery error. The class
    % remains a full subclass of image_vector; only these specific methods are
    % masked. See docs/fmri_surface_data_methods.md ("Unsupported methods").
    % ---------------------------------------------------------------------
    methods (Hidden)
        function varargout = flip(varargin),                        unsupported_surface_method('flip'); varargout = {}; end
        function varargout = isosurface(varargin),                  unsupported_surface_method('isosurface'); varargout = {}; end
        function varargout = interpolate(varargin),                 unsupported_surface_method('interpolate'); varargout = {}; end
        function varargout = resample_space(varargin),              unsupported_surface_method('resample_space', 'resample_space (surface): use resample_space on to_fmri_data(obj) for the subcortex, or keep surface data in its native grayordinate space'); varargout = {}; end
        function varargout = resample_space_simple_reference(varargin), unsupported_surface_method('resample_space_simple_reference'); varargout = {}; end
        function varargout = resample_time(varargin),               unsupported_surface_method('resample_time'); varargout = {}; end
        function varargout = extract_gray_white_csf(varargin),      unsupported_surface_method('extract_gray_white_csf'); varargout = {}; end
        function varargout = slice_movie(varargin),                 unsupported_surface_method('slice_movie'); varargout = {}; end
        function varargout = display_slices(varargin),              unsupported_surface_method('display_slices'); varargout = {}; end
        function varargout = rmssd_movie(varargin),                 unsupported_surface_method('rmssd_movie'); varargout = {}; end
        function varargout = render_on_cerebellar_flatmap(varargin),unsupported_surface_method('render_on_cerebellar_flatmap'); varargout = {}; end
        function varargout = rebuild_volinfo_from_dat(varargin),    unsupported_surface_method('rebuild_volinfo_from_dat'); varargout = {}; end
        function varargout = trim_mask(varargin),                   unsupported_surface_method('trim_mask'); varargout = {}; end
        function varargout = define_space_mapping(varargin),        unsupported_surface_method('define_space_mapping'); varargout = {}; end
        function varargout = plot_current_orthviews_coord(varargin),unsupported_surface_method('plot_current_orthviews_coord'); varargout = {}; end
        function varargout = read_from_file(varargin),              unsupported_surface_method('read_from_file', 'read_from_file: construct with fmri_surface_data(filename) instead (CIFTI/GIFTI are read natively)'); varargout = {}; end
        function varargout = check_image_filenames(varargin),       unsupported_surface_method('check_image_filenames'); varargout = {}; end
        function varargout = searchlight(varargin),                 unsupported_surface_method('searchlight'); varargout = {}; end
        function varargout = extract_roi_averages(varargin),        unsupported_surface_method('extract_roi_averages', 'extract_roi_averages (volume ROI): use apply_parcellation(obj, surface_atlas) for surface data'); varargout = {}; end
        function varargout = subdivide_by_atlas(varargin),          unsupported_surface_method('subdivide_by_atlas'); varargout = {}; end
        function varargout = expand_into_atlas_subregions(varargin),unsupported_surface_method('expand_into_atlas_subregions'); varargout = {}; end
        function varargout = table_of_atlas_regions_covered(varargin), unsupported_surface_method('table_of_atlas_regions_covered'); varargout = {}; end
        function varargout = pattern_surf_plot_mip(varargin),       unsupported_surface_method('pattern_surf_plot_mip'); varargout = {}; end
    end

end % classdef


% =========================================================================
function unsupported_surface_method(name, custom_msg)
% Consistent, informative error for volume-only methods masked on surface data.
if nargin >= 2 && ~isempty(custom_msg)
    hint = custom_msg;
else
    hint = ['For cortical surface data use surface(obj) / render_on_surface; ' ...
            'for the subcortical volume use to_fmri_data(obj) and the fmri_data method.'];
end
error('fmri_surface_data:unsupportedMethod', ...
    ['%s is a volume-only image_vector method and is not supported for ' ...
     'fmri_surface_data (surface/grayordinate) objects.\n%s'], name, hint);
end


% =========================================================================
% Local helper functions (private to the constructor)
% =========================================================================

function obj = load_surface_file(obj, fname)
% Dispatch a file to the native CIFTI or GIFTI reader based on extension.
if exist(fname, 'file') ~= 2
    error('fmri_surface_data:notfound', 'File not found: %s', fname);
end
lc = lower(fname);
if endsWith(lc, '.gii')
    g = canlab_read_gifti(fname);
    obj = from_gifti_struct(obj, g);
elseif endsWith(lc, '.nii') || endsWith(lc, '.nii.gz')
    cii = canlab_read_cifti(fname);
    obj = from_cifti_struct(obj, cii);
else
    error('fmri_surface_data:badext', ...
        'Unrecognized surface/grayordinate extension: %s (expected .nii CIFTI or .gii GIFTI).', fname);
end
obj.fullpath = fname;
end


% -------------------------------------------------------------------------
function obj = from_cifti_struct(obj, cii)
% Build the object from a canlab_read_cifti output struct.
obj.dat = single(cii.cdata);
obj.imagetype = cii.intent;

bm = cii.diminfo{1};                          % dense dimension
if ~isfield(bm, 'cluster'), bm.cluster = []; end
bm.grayordinate_type = infer_grayord_type(bm);
obj.brain_model = bm;
obj.surface_space = infer_surface_space(bm);

% Populate the inherited volInfo slot for the subcortical voxel sub-block only
% (empty for surface-only objects). brain_model remains the geometry truth.
obj.volInfo = build_volinfo_subblock(bm);

% Maps (second) dimension -> per-map annotations
md = cii.diminfo{2};
switch md.type
    case 'scalars'
        obj.image_names = reshape({md.maps.name}, [], 1);
    case 'labels'
        obj.image_names = reshape({md.maps.name}, [], 1);
        obj.label_table = struct2labeltable(md.maps(1).table);
        obj.additional_info.label_tables = {md.maps.table};
    case 'series'
        obj.series_info = struct('start', md.seriesStart, 'step', md.seriesStep, ...
            'unit', md.seriesUnit, 'exponent', getfielddef(md, 'seriesExponent', 0));
        obj.image_names = arrayfun(@(k) sprintf('t%d', k), (1:size(obj.dat,2))', 'unif', 0);
end

obj.removed_voxels = false(size(obj.dat, 1), 1);
obj.removed_images = false(size(obj.dat, 2), 1);

% Stash source XML/header so write() can faithfully re-emit if unchanged
if isfield(cii, 'xml'), obj.additional_info.cifti_xml = cii.xml; end
if isfield(cii, 'hdr'), obj.additional_info.cifti_hdr = cii.hdr; end

obj.history{end+1} = sprintf(['fmri_surface_data created from CIFTI (%s): ' ...
    '%d grayordinates x %d maps, space=%s'], cii.intent, size(obj.dat,1), size(obj.dat,2), obj.surface_space);
end


% -------------------------------------------------------------------------
function obj = from_gifti_struct(obj, g)
% Build from a canlab_read_gifti output struct (functional/label or geometry).
hasData = isfield(g, 'cdata') && ~isempty(g.cdata) && ~iscell(g.cdata);
hasGeom = isfield(g, 'vertices') && ~isempty(g.vertices);

if hasData
    obj.dat = single(g.cdata);
    n = size(g.cdata, 1);
    m = struct('struct', 'CORTEX', 'type', 'surf', 'start', 1, 'count', n, ...
        'numvert', n, 'vertlist', 0:n-1, 'voxlist', []);
    bm = struct('type', 'dense', 'length', n, 'models', {{m}}, 'vol', []);
    bm.grayordinate_type = 'cortex_only';
    bm.cluster = [];
    obj.brain_model = bm;
    obj.surface_space = infer_surface_space(bm);
    obj.imagetype = 'func';
    if isfield(g, 'labels') && ~isempty(g.labels)
        obj.label_table = struct2labeltable(g.labels);
        obj.imagetype = 'label';
    end
    if isfield(g, 'intents') && ~isempty(g.intents)
        obj.image_names = reshape(g.intents, [], 1);
    end
    obj.removed_voxels = false(n, 1);
    obj.removed_images = false(size(obj.dat, 2), 1);
    obj.history{end+1} = sprintf('fmri_surface_data created from GIFTI (%s): %d vertices x %d maps', ...
        obj.imagetype, n, size(obj.dat,2));
end

if hasGeom
    % Store geometry; if this is a pure .surf.gii (no data) the object holds a mesh.
    obj.geom = struct('vertices', g.vertices, 'faces', g.faces);
    if ~hasData
        obj.history{end+1} = sprintf('fmri_surface_data loaded surface geometry: %d vertices, %d faces', ...
            size(g.vertices,1), size(g.faces,1));
    end
end
end


% -------------------------------------------------------------------------
function obj = copy_matching_fields(obj, s)
p = properties(obj);
fn = fieldnames(s);
for i = 1:numel(fn)
    if ismember(fn{i}, p)
        obj.(fn{i}) = s.(fn{i});
    end
end
if ~isempty(obj.dat), obj.dat = single(obj.dat); end
obj = normalize_removed(obj);
end


% -------------------------------------------------------------------------
function obj = normalize_removed(obj)
% Keep the vestigial removed_* vectors all-false and length-correct (D5b),
% whenever .dat is set via key-value / struct / recast construction.
if ~isempty(obj.dat)
    if numel(obj.removed_voxels) ~= size(obj.dat, 1)
        obj.removed_voxels = false(size(obj.dat, 1), 1);
    end
    if numel(obj.removed_images) ~= size(obj.dat, 2)
        obj.removed_images = false(size(obj.dat, 2), 1);
    end
end
end


% -------------------------------------------------------------------------
function obj = recast_from_object(obj, other)
p = properties(obj);
op = properties(other);
for i = 1:numel(op)
    if ismember(op{i}, p)
        try, obj.(op{i}) = other.(op{i}); catch, end %#ok<NOCOM>
    end
end
if ~isempty(obj.dat), obj.dat = single(obj.dat); end
obj = normalize_removed(obj);
obj.history{end+1} = sprintf('Recast from %s', class(other));
end


% -------------------------------------------------------------------------
function obj = parse_keyvalue(obj, args)
p = properties(obj);
i = 1;
while i <= numel(args)
    key = args{i};
    if (ischar(key) || isstring(key)) && ismember(char(key), p)
        obj.(char(key)) = args{i+1};
        i = i + 2;
    elseif (ischar(key) || isstring(key)) && (contains(char(key), '.nii') || contains(char(key), '.gii')) ...
            && exist(char(key), 'file')
        obj.fullpath = char(key);
        i = i + 1;
    else
        warning('fmri_surface_data:badinput', 'Unknown field or argument: %s', char(string(key)));
        i = i + 1;
    end
end
if ~isempty(obj.dat), obj.dat = single(obj.dat); end
obj = normalize_removed(obj);
end


% -------------------------------------------------------------------------
function s = infer_surface_space(bm)
% Infer a canonical surface-space tag from the per-hemisphere vertex count.
s = '';
nv = NaN;
for i = 1:numel(bm.models)
    if strcmp(bm.models{i}.type, 'surf') && ~isempty(bm.models{i}.numvert) ...
            && ~isnan(bm.models{i}.numvert)
        nv = bm.models{i}.numvert; break
    end
end
switch nv
    case 32492,  s = 'fsLR_32k';
    case 163842, s = 'fsaverage_164k';
    case 40962,  s = 'fsaverage6';
    case 10242,  s = 'fsaverage5';
    case 2562,   s = 'fsaverage4';
    otherwise
        if ~isnan(nv), s = sprintf('surf_%dverts', nv); end
end
end


% -------------------------------------------------------------------------
function t = infer_grayord_type(bm)
% '91k'-style grayordinate (surface + subcortical voxels) vs cortex-only.
hasvox = any(cellfun(@(m) strcmp(m.type, 'vox'), bm.models));
hassurf = any(cellfun(@(m) strcmp(m.type, 'surf'), bm.models));
if hasvox && hassurf
    if bm.length == 91282, t = '91k';
    else, t = sprintf('grayordinate_%d', bm.length);
    end
elseif hassurf
    t = 'cortex_only';
elseif hasvox
    t = 'volume_only';
else
    t = 'unknown';
end
end


% -------------------------------------------------------------------------
function v = getfielddef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)) && ~(isnumeric(s.(f)) && all(isnan(s.(f)(:))))
    v = s.(f);
else
    v = d;
end
end
