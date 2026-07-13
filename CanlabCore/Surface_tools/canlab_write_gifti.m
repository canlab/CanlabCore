function canlab_write_gifti(filename, gii, varargin)
% Write a GIFTI (.gii) file natively in MATLAB, with no external toolbox.
%
% :Usage:
% ::
%     canlab_write_gifti(filename, gii)
%     canlab_write_gifti(filename, gii, 'encoding', 'GZipBase64Binary')
%
% Writes geometry (vertices+faces), per-vertex functional/shape data, and/or a
% label table to a GIFTI file. The struct layout matches canlab_read_gifti's
% output, so read->write->read round-trips exactly. Depends only on core MATLAB
% plus the JVM (for zlib deflate) -- no SPM/FieldTrip gifti toolbox required.
%
% :Inputs:
%
%   **filename:**
%        Char/string. Output path (e.g. '/tmp/out.surf.gii').
%
%   **gii:**
%        Struct (as from canlab_read_gifti) with any of:
%          .vertices : [N x 3] vertex coords -> NIFTI_INTENT_POINTSET
%          .faces    : [M x 3] 1-based faces -> NIFTI_INTENT_TRIANGLE (written 0-based)
%          .cdata    : [N x K] per-vertex data -> one DataArray per column
%          .intents  : 1xK cell of intent strings for cdata columns
%                      (default NIFTI_INTENT_NONE, or _LABEL if .labels present)
%          .labels   : struct array (.key,.name,.rgba) -> <LabelTable>
%
% :Optional Inputs:
%   **'encoding':** 'GZipBase64Binary' (default) | 'Base64Binary' | 'ASCII'
%
% :Examples:
% ::
%     g = canlab_read_gifti(which('S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii'));
%     canlab_write_gifti('/tmp/copy.surf.gii', g);
%
% :See also: canlab_read_gifti, canlab_read_cifti, canlab_write_cifti

encoding = 'GZipBase64Binary';
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'encoding', encoding = varargin{i+1};
        otherwise, error('canlab_write_gifti:badopt', 'Unknown option: %s', varargin{i});
    end
end

if ~isstruct(gii)
    error('canlab_write_gifti:input', 'Second argument must be a struct (see help).');
end
gii = local_fill_defaults(gii);

% Count data arrays
nda = 0;
if ~isempty(gii.vertices), nda = nda + 1; end
if ~isempty(gii.faces),    nda = nda + 1; end
if ~isempty(gii.cdata),    nda = nda + size(gii.cdata, 2); end
if nda == 0
    error('canlab_write_gifti:empty', 'Nothing to write (no vertices/faces/cdata).');
end

fid = fopen(char(filename), 'w');
if fid < 0, error('canlab_write_gifti:open', 'Cannot open for writing: %s', filename); end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid, '<!DOCTYPE GIFTI SYSTEM "http://www.nitrc.org/frs/download.php/1594/gifti.dtd">\n');
fprintf(fid, '<GIFTI Version="1.0" NumberOfDataArrays="%d">\n', nda);

% Optional label table
if ~isempty(gii.labels)
    local_write_labeltable(fid, gii.labels);
end

% Geometry
if ~isempty(gii.vertices)
    local_write_da(fid, 'NIFTI_INTENT_POINTSET', 'NIFTI_TYPE_FLOAT32', ...
        single(gii.vertices), encoding);
end
if ~isempty(gii.faces)
    local_write_da(fid, 'NIFTI_INTENT_TRIANGLE', 'NIFTI_TYPE_INT32', ...
        int32(gii.faces - 1), encoding);   % write 0-based
end

% Per-vertex data: one DataArray per column
for c = 1:size(gii.cdata, 2)
    intent = 'NIFTI_INTENT_NONE';
    if c <= numel(gii.intents) && ~isempty(gii.intents{c})
        intent = gii.intents{c};
    end
    col = gii.cdata(:, c);
    if contains(intent, 'LABEL')
        local_write_da(fid, intent, 'NIFTI_TYPE_INT32', int32(col), encoding);
    else
        local_write_da(fid, intent, 'NIFTI_TYPE_FLOAT32', single(col), encoding);
    end
end

fprintf(fid, '</GIFTI>\n');
end % main


% -------------------------------------------------------------------------
function gii = local_fill_defaults(gii)
flds = {'vertices', 'faces', 'cdata', 'intents', 'labels'};
for i = 1:numel(flds)
    if ~isfield(gii, flds{i}), gii.(flds{i}) = []; end
end
if isempty(gii.intents), gii.intents = {}; end
if ~isempty(gii.labels) && ~isempty(gii.cdata) && isempty(gii.intents)
    gii.intents = repmat({'NIFTI_INTENT_LABEL'}, 1, size(gii.cdata, 2));
end
end


% -------------------------------------------------------------------------
function local_write_da(fid, intent, dtype, data, encoding)
% data is [Dim0 x Dim1] in MATLAB column-major; GIFTI stores RowMajorOrder.
d0 = size(data, 1);
d1 = size(data, 2);
dimn = 1 + (d1 > 1);

% Row-major byte stream = bytes of the transpose, read column-major
dt = data';
switch dtype
    case 'NIFTI_TYPE_FLOAT32', bytes = typecast(single(dt(:))', 'uint8');
    case 'NIFTI_TYPE_INT32',   bytes = typecast(int32(dt(:))',  'uint8');
    case 'NIFTI_TYPE_FLOAT64', bytes = typecast(double(dt(:))', 'uint8');
    case 'NIFTI_TYPE_UINT8',   bytes = uint8(dt(:))';
    otherwise, error('canlab_write_gifti:dtype', 'Unsupported DataType %s', dtype);
end

if dimn == 1
    dimattr = sprintf('Dimensionality="1" Dim0="%d"', d0);
else
    dimattr = sprintf('Dimensionality="2" Dim0="%d" Dim1="%d"', d0, d1);
end

fprintf(fid, ['<DataArray Intent="%s" DataType="%s" ArrayIndexingOrder="RowMajorOrder" ' ...
    '%s Encoding="%s" Endian="LittleEndian" ExternalFileName="" ExternalFileOffset="">\n'], ...
    intent, dtype, dimattr, encoding);

switch encoding
    case 'ASCII'
        if contains(dtype, 'INT')
            fmt = '%d ';
        else
            fmt = '%.9g ';   % full single-precision width
        end
        fprintf(fid, '<Data>');
        fprintf(fid, fmt, double(dt(:)));
        fprintf(fid, '</Data>\n');
    case 'Base64Binary'
        fprintf(fid, '<Data>%s</Data>\n', char(matlab.net.base64encode(bytes)));
    case 'GZipBase64Binary'
        comp = canlab_zlib_deflate(bytes);
        fprintf(fid, '<Data>%s</Data>\n', char(matlab.net.base64encode(comp)));
    otherwise
        error('canlab_write_gifti:encoding', 'Unknown encoding %s', encoding);
end

fprintf(fid, '</DataArray>\n');
end


% -------------------------------------------------------------------------
function local_write_labeltable(fid, labels)
fprintf(fid, '<LabelTable>\n');
for i = 1:numel(labels)
    rgba = labels(i).rgba;
    if numel(rgba) < 4, rgba = [rgba(:)' ones(1, 4 - numel(rgba))]; end
    nm = labels(i).name;
    if isempty(nm), nm = sprintf('label_%d', labels(i).key); end
    fprintf(fid, '<Label Key="%d" Red="%g" Green="%g" Blue="%g" Alpha="%g"><![CDATA[%s]]></Label>\n', ...
        labels(i).key, rgba(1), rgba(2), rgba(3), rgba(4), nm);
end
fprintf(fid, '</LabelTable>\n');
end
