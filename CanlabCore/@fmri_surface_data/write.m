function write(obj, varargin)
% write Write an fmri_surface_data object to a CIFTI-2 (.nii) or GIFTI (.gii) file.
%
% :Usage:
% ::
%     write(obj)                          % writes to obj.fullpath
%     write(obj, '/path/out.dscalar.nii') % CIFTI-2 (dscalar/dtseries/dlabel)
%     write(obj, 'fname', '/path/out.func.gii')   % GIFTI
%
% Dispatches on the output extension to the native writers canlab_write_cifti /
% canlab_write_gifti (M1) -- no external toolbox is required. The CIFTI maps
% dimension is rebuilt from the object's intent + image_names / label_table /
% series_info. If the object was loaded from a CIFTI and its grayordinate layout
% is unchanged, the original CIFTI XML is re-emitted faithfully; otherwise the
% XML is regenerated from brain_model.
%
% :Inputs:
%   **obj:**       an fmri_surface_data object.
%   **filename:**  output path (positional) or via 'fname'/'filename' keyword.
%                  Defaults to obj.fullpath.
%
% :See also: canlab_write_cifti, canlab_write_gifti, fmri_surface_data, to_fmri_data

% ---- Resolve filename ----
fname = '';
i = 1;
while i <= numel(varargin)
    a = varargin{i};
    if (ischar(a) || isstring(a)) && any(strcmpi(char(a), {'fname','filename'})) && i < numel(varargin)
        fname = char(varargin{i+1}); i = i + 2;
    elseif (ischar(a) || isstring(a)) && (contains(char(a), '.nii') || contains(char(a), '.gii'))
        fname = char(a); i = i + 1;
    else
        i = i + 1;
    end
end
if isempty(fname), fname = char(obj.fullpath); end
if isempty(fname)
    error('fmri_surface_data:write:nofname', ...
        'No output filename. Pass one or set obj.fullpath.');
end

lc = lower(fname);

% ===================== GIFTI =====================
if endsWith(lc, '.gii')
    g = struct('vertices', [], 'faces', [], 'cdata', [], 'intents', {{}}, 'labels', []);
    if ~isempty(obj.geom) && isfield(obj.geom, 'vertices')
        g.vertices = obj.geom.vertices;
        g.faces = obj.geom.faces;
    end
    if ~isempty(obj.dat)
        g.cdata = double(obj.dat);
        islabel = contains(lc, '.label.') || strcmp(obj.intent, 'label') || strcmp(obj.intent, 'dlabel');
        if islabel
            g.intents = repmat({'NIFTI_INTENT_LABEL'}, 1, size(g.cdata, 2));
            g.labels = obj.label_table;
        else
            g.intents = repmat({'NIFTI_INTENT_NONE'}, 1, size(g.cdata, 2));
        end
    end
    canlab_write_gifti(fname, g);
    return
end

% ===================== CIFTI-2 =====================
intent = local_resolve_intent(obj.intent, lc);

cii = struct();
cii.cdata = double(obj.dat);
cii.intent = intent;
cii.diminfo = {obj.brain_model, local_build_maps_dim(obj, intent)};

% Faithful re-emit if the stashed source XML still matches the current layout
ai = obj.additional_info;
if isstruct(ai) && isfield(ai, 'cifti_xml') && isfield(ai, 'cifti_hdr') ...
        && ~isempty(ai.cifti_xml) && ~isempty(ai.cifti_hdr)
    h = ai.cifti_hdr;
    if isfield(h, 'dim') && numel(h.dim) >= 7
        gray = size(cii.cdata, 1); maps = size(cii.cdata, 2);
        if (h.dim(6) == maps && h.dim(7) == gray) || (h.dim(6) == gray && h.dim(7) == maps)
            cii.xml = ai.cifti_xml;
            cii.hdr = h;
        end
    end
end

canlab_write_cifti(fname, cii);
end


% =========================================================================
function intent = local_resolve_intent(objintent, lc)
if ismember(objintent, {'dscalar','dtseries','dlabel','dconn'})
    intent = objintent;
elseif contains(lc, '.dscalar.'), intent = 'dscalar';
elseif contains(lc, '.dtseries.'), intent = 'dtseries';
elseif contains(lc, '.dlabel.'), intent = 'dlabel';
elseif strcmp(objintent, 'label'), intent = 'dlabel';
elseif strcmp(objintent, 'func') || strcmp(objintent, 'shape'), intent = 'dscalar';
else, intent = 'dscalar';
end
end


% =========================================================================
function md = local_build_maps_dim(obj, intent)
nMaps = size(obj.dat, 2);

names = obj.image_names;
if isempty(names) || numel(names) ~= nMaps
    names = arrayfun(@(k) sprintf('map_%d', k), 1:nMaps, 'UniformOutput', false);
end
names = reshape(names, 1, []);

switch intent
    case 'dlabel'
        % Per-map label tables: prefer stashed per-map tables, else reuse one
        tables = {};
        if isstruct(obj.additional_info) && isfield(obj.additional_info, 'label_tables')
            tables = obj.additional_info.label_tables;
        end
        maps = struct('name', {}, 'table', {});
        for k = 1:nMaps
            if k <= numel(tables) && ~isempty(tables{k})
                tbl = tables{k};
            else
                tbl = obj.label_table;
            end
            maps(k) = struct('name', names{k}, 'table', tbl);
        end
        md = struct('type', 'labels', 'length', nMaps, 'maps', maps);

    case 'dtseries'
        si = obj.series_info;
        md = struct('type', 'series', 'length', nMaps);
        md.seriesStart    = local_getdef(si, 'start', 0);
        md.seriesStep     = local_getdef(si, 'step', 1);
        md.seriesUnit     = local_getstr(si, 'unit', 'SECOND');
        md.seriesExponent = local_getdef(si, 'exponent', 0);

    otherwise   % dscalar (and fallback)
        maps = struct('name', {}, 'table', {});
        for k = 1:nMaps
            maps(k) = struct('name', names{k}, 'table', []);
        end
        md = struct('type', 'scalars', 'length', nMaps, 'maps', maps);
end
end


function v = local_getdef(s, f, d)
v = d;
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    val = double(s.(f));
    if ~all(isnan(val(:)))
        v = s.(f);
    end
end
end

function v = local_getstr(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end
