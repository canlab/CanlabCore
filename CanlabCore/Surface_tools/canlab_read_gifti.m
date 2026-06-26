function gii = canlab_read_gifti(filename)
% Read a GIFTI (.gii) file natively in MATLAB, with no external toolbox.
%
% :Usage:
% ::
%     gii = canlab_read_gifti(filename)
%
% Reads geometry (.surf.gii), functional/shape (.func.gii / .shape.gii) and
% label (.label.gii) GIFTI files. GIFTI is a UTF-8 XML container holding one or
% more <DataArray> elements; payloads may be ASCII, Base64Binary, or
% GZipBase64Binary (base64 of a raw zlib stream). This function depends only on
% core MATLAB plus the JVM (for zlib inflate of GZipBase64Binary data) -- it does
% NOT require the SPM/FieldTrip gifti toolbox or Connectome Workbench.
%
% :Inputs:
%
%   **filename:**
%        Char/string. Full path to a .gii file (.surf/.func/.shape/.label.gii).
%
% :Outputs:
%
%   **gii:**
%        Struct with fields:
%          .vertices : [N x 3] double vertex mm coordinates (POINTSET), or []
%          .faces    : [M x 3] double, 1-based triangle indices (TRIANGLE), or []
%          .cdata    : [N x K] per-vertex data, K columns from K data arrays, or []
%          .intents  : 1xK cell of NIFTI_INTENT_* strings for the cdata columns
%          .labels   : struct array (.key,.name,.rgba) from a GIFTI LabelTable, or []
%          .gifti_version : char, from the root <GIFTI Version="..."> attribute
%
% :Examples:
% ::
%     g = canlab_read_gifti(which('S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii'));
%     patch('Faces', g.faces, 'Vertices', g.vertices, 'FaceColor', [.7 .7 .7]); axis equal
%
%     lbl = canlab_read_gifti('Glasser_2016.32k.L.label.gii');
%     unique(lbl.cdata)              % integer region keys
%     {lbl.labels.name}             % region names
%
% :See also:
%   canlab_write_gifti, canlab_read_cifti, canlab_write_cifti, fmri_surface_data

% ..
%    Part of the CANlab fmri_surface_data surface-data toolset. See
%    docs/fmri_surface_data_design_plan.md. No external toolbox required.
% ..

if nargin < 1 || isempty(filename)
    error('canlab_read_gifti:input', 'Provide a path to a .gii file.');
end
filename = char(filename);
if exist(filename, 'file') ~= 2
    error('canlab_read_gifti:notfound', 'File not found: %s', filename);
end

txt = fileread(filename);

gii = struct('vertices', [], 'faces', [], 'cdata', [], 'intents', {{}}, ...
             'labels', [], 'gifti_version', '');

% Root version
v = regexp(txt, '<GIFTI[^>]*\sVersion="([^"]*)"', 'tokens', 'once');
if ~isempty(v), gii.gifti_version = v{1}; end

% Label table (optional; .label.gii)
gii.labels = local_parse_labeltable(txt);

% Each DataArray
das = regexp(txt, '(?s)<DataArray\s.*?</DataArray>', 'match');
cdata = {};
for k = 1:numel(das)
    [M, intent] = local_parse_dataarray(das{k});
    if isempty(M), continue; end
    if contains(intent, 'POINTSET')
        gii.vertices = M;
    elseif contains(intent, 'TRIANGLE')
        gii.faces = M + 1;            % store 1-based for patch()
    else
        cdata{end+1} = M(:);          %#ok<AGROW> per-vertex data column
        gii.intents{end+1} = intent;  %#ok<AGROW>
    end
end

if ~isempty(cdata)
    n = cellfun(@numel, cdata);
    if numel(unique(n)) == 1
        gii.cdata = double([cdata{:}]);
    else
        % Ragged data arrays: keep as cell to avoid silent truncation
        gii.cdata = cdata;
        warning('canlab_read_gifti:ragged', ...
            'DataArrays have differing lengths; returning .cdata as a cell array.');
    end
end

end % main function


% -------------------------------------------------------------------------
function [M, intent] = local_parse_dataarray(da)
% Parse one <DataArray>...</DataArray> block into a numeric matrix.

ga = @(name) local_attr(da, name);
intent = ga('Intent');
dtype  = ga('DataType');
order  = ga('ArrayIndexingOrder');
enc    = ga('Encoding');
endian = ga('Endian');

dimn = str2double(ga('Dimensionality'));
if isnan(dimn) || dimn < 1, dimn = 1; end
dims = ones(1, dimn);
for d = 0:dimn-1
    dims(d+1) = str2double(ga(sprintf('Dim%d', d)));
end

switch dtype
    case 'NIFTI_TYPE_FLOAT32', cls = 'single';
    case 'NIFTI_TYPE_FLOAT64', cls = 'double';
    case 'NIFTI_TYPE_INT32',   cls = 'int32';
    case 'NIFTI_TYPE_UINT32',  cls = 'uint32';
    case 'NIFTI_TYPE_INT16',   cls = 'int16';
    case 'NIFTI_TYPE_UINT8',   cls = 'uint8';
    case 'NIFTI_TYPE_INT8',    cls = 'int8';
    otherwise,                 cls = 'single';
end

payload = regexp(da, '(?s)<Data>(.*?)</Data>', 'tokens', 'once');
if isempty(payload), M = []; return; end
data_str = payload{1};

if strcmp(enc, 'ASCII')
    vals = cast(sscanf(data_str, '%f'), cls);   % whitespace-delimited numbers
else
    b64 = data_str(~isspace(data_str));         % base64 must have whitespace stripped
    if isempty(b64), M = []; return; end
    raw = matlab.net.base64decode(b64);
    if strcmp(enc, 'GZipBase64Binary')
        raw = canlab_zlib_inflate(raw);
    end
    vals = typecast(raw(:), cls);
    if strcmpi(endian, 'BigEndian'), vals = swapbytes(vals); end
end
vals = double(vals);

if numel(dims) >= 2 && dims(2) > 1
    if strcmp(order, 'RowMajorOrder')
        M = reshape(vals, [dims(2) dims(1)])';
    else
        M = reshape(vals, dims);
    end
else
    M = vals(:);
end
end


% -------------------------------------------------------------------------
function labels = local_parse_labeltable(txt)
% Parse an optional GIFTI <LabelTable> into a struct array (key,name,rgba).

labels = [];
lt = regexp(txt, '(?s)<LabelTable>(.*?)</LabelTable>', 'tokens', 'once');
if isempty(lt), return; end
items = regexp(lt{1}, '(?s)<Label\s.*?</Label>', 'match');
if isempty(items), return; end

labels = struct('key', {}, 'name', {}, 'rgba', {});
for i = 1:numel(items)
    it = items{i};
    key = str2double(local_attr(it, 'Key'));
    r = str2double(local_attr(it, 'Red'));
    g = str2double(local_attr(it, 'Green'));
    b = str2double(local_attr(it, 'Blue'));
    a = str2double(local_attr(it, 'Alpha'));
    nm = regexp(it, '(?s)<Label\s[^>]*>(.*?)</Label>', 'tokens', 'once');
    name = '';
    if ~isempty(nm)
        name = strtrim(nm{1});
        name = regexprep(name, '<!\[CDATA\[(.*?)\]\]>', '$1');
    end
    labels(end+1) = struct('key', key, 'name', name, ...
        'rgba', [r g b a]); %#ok<AGROW>
end
end


% -------------------------------------------------------------------------
function val = local_attr(str, name)
% Extract a single XML attribute value (or '' if absent).
tok = regexp(str, [name '\s*=\s*"([^"]*)"'], 'tokens', 'once');
if isempty(tok), val = ''; else, val = tok{1}; end
end
