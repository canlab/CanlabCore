function cii = canlab_read_cifti(filename)
% Read a CIFTI-2 file (.dscalar/.dtseries/.dlabel/.dconn.nii) natively in MATLAB.
%
% :Usage:
% ::
%     cii = canlab_read_cifti(filename)
%
% A CIFTI-2 file is a NIfTI-2 container (540-byte little-endian header) whose
% header extension with ecode==32 holds an XML description of the matrix
% (BrainModels = surface vertices + subcortical voxels; plus a Scalars/Series/
% Labels map for the second dimension). The data matrix follows at vox_offset.
% This reader uses only fopen/fread + regexp XML parsing -- it does NOT require
% cifti-matlab, FieldTrip, gifti, or Connectome Workbench at runtime.
%
% Only modern CIFTI-2 is supported (legacy CIFTI-1 requires wb_command to
% upgrade). The output field names mirror cifti-matlab's so existing helpers
% such as get_cifti_data.m can consume it unchanged.
%
% :Inputs:
%
%   **filename:**
%        Char/string. Path to a .dscalar.nii / .dtseries.nii / .dlabel.nii /
%        .dconn.nii (CIFTI-2) file.
%
% :Outputs:
%
%   **cii:**
%        Struct with fields:
%          .cdata    : [nGrayordinates x nMaps] numeric data matrix
%          .intent   : 'dscalar' | 'dtseries' | 'dlabel' | 'dconn' | 'unknown'
%          .diminfo  : 1x2 cell. diminfo{1} describes the grayordinate (dense)
%                      dimension; diminfo{2} the map/series/label dimension.
%             diminfo{1}.type = 'dense'
%             diminfo{1}.length
%             diminfo{1}.models{i}: .struct (e.g. 'CORTEX_LEFT'), .type 'surf'|'vox',
%                  .start (1-based row), .count, .numvert (surf), .vertlist (0-based),
%                  .voxlist (3xN IJK, 0-based)
%             diminfo{1}.vol.dims (1x3), diminfo{1}.vol.sform (4x4, mm)
%             diminfo{2}.type = 'scalars'|'series'|'labels'
%             diminfo{2}.maps(k).name (+ .table for labels: .key,.name,.rgba)
%             diminfo{2}.seriesStart/.seriesStep/.seriesUnit (series)
%          .hdr      : parsed NIfTI-2 header fields (dim, datatype, vox_offset, ...)
%          .xml      : raw CIFTI XML char array (used by canlab_write_cifti)
%
% :Examples:
% ::
%     c = canlab_read_cifti(which('transcriptomic_gradients.dscalar.nii'));
%     size(c.cdata)                         % [grayordinates x maps]
%     cellfun(@(m) m.struct, c.diminfo{1}.models, 'unif', 0)
%
% :See also:
%   canlab_write_cifti, canlab_read_gifti, get_cifti_data, fmri_surface_data

% ..
%    Part of the CANlab fmri_surface_data surface-data toolset. See
%    docs/fmri_surface_data_design_plan.md. No external toolbox required.
% ..

if nargin < 1 || isempty(filename)
    error('canlab_read_cifti:input', 'Provide a path to a CIFTI-2 .nii file.');
end
filename = char(filename);
if exist(filename, 'file') ~= 2
    error('canlab_read_cifti:notfound', 'File not found: %s', filename);
end

% ---- Open and detect endianness via sizeof_hdr (must be 540 for NIfTI-2) ----
fid = fopen(filename, 'r', 'l');
if fid < 0, error('canlab_read_cifti:open', 'Cannot open %s', filename); end
cleanupObj = onCleanup(@() fclose(fid));

sizeof_hdr = fread(fid, 1, 'int32');
machine = 'l';
if sizeof_hdr ~= 540
    fclose(fid);
    fid = fopen(filename, 'r', 'b'); machine = 'b';
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    sizeof_hdr = fread(fid, 1, 'int32');
    if sizeof_hdr ~= 540
        error('canlab_read_cifti:notnifti2', ...
            ['sizeof_hdr=%d (expected 540). Not a NIfTI-2/CIFTI-2 file. ' ...
             'Legacy CIFTI-1 is not supported.'], sizeof_hdr);
    end
end

% ---- NIfTI-2 header (fixed offsets) ----
hdr.magic      = fread(fid, 8, '*char')';      % @4
hdr.datatype   = fread(fid, 1, 'int16');       % @12
hdr.bitpix     = fread(fid, 1, 'int16');       % @14
hdr.dim        = fread(fid, 8, 'int64')';      % @16
hdr.intent_p   = fread(fid, 3, 'double')';     % @80
hdr.pixdim     = fread(fid, 8, 'double')';     % @104
hdr.vox_offset = fread(fid, 1, 'int64');       % @168
hdr.scl_slope  = fread(fid, 1, 'double');      % @176
hdr.scl_inter  = fread(fid, 1, 'double');      % @184
fseek(fid, 504, 'bof');
hdr.intent_code = fread(fid, 1, 'int32');      % @504
hdr.intent_name = deblank(fread(fid, 16, '*char')'); % @508
hdr.machine = machine;

% ---- Header extension: find the ecode==32 CIFTI XML ----
fseek(fid, 540, 'bof');
extender = fread(fid, 4, 'uint8');
xml = '';
if ~isempty(extender) && extender(1) ~= 0
    while ftell(fid) < hdr.vox_offset
        esize = fread(fid, 1, 'int32');
        ecode = fread(fid, 1, 'int32');
        if isempty(esize) || esize <= 8, break; end
        edata = fread(fid, esize - 8, '*char')';
        if ecode == 32, xml = edata; end
    end
end
if isempty(xml)
    error('canlab_read_cifti:noxml', 'No CIFTI XML extension (ecode 32) found.');
end
% Trim trailing NULs/whitespace padding from the XML blob
xml = regexprep(xml, '\x00+$', '');
cii.xml = xml;

% ---- Read the data matrix ----
rowLen = hdr.dim(6);    % NIfTI dim[5]
nRows  = hdr.dim(7);    % NIfTI dim[6]
prec = local_nifti_precision(hdr.datatype);
fseek(fid, hdr.vox_offset, 'bof');
v = fread(fid, rowLen * nRows, [prec '=>' prec]);
v = double(v);
if isfinite(hdr.scl_slope) && hdr.scl_slope ~= 0 && ...
        ~(hdr.scl_slope == 1 && hdr.scl_inter == 0)
    v = v * hdr.scl_slope + hdr.scl_inter;
end

% ---- Parse XML into diminfo ----
[diminfo, intent] = local_parse_cifti_xml(xml, hdr.intent_code, hdr.intent_name);

% ---- Orient cdata to [grayordinates x maps] robustly ----
G = 0;
if ~isempty(diminfo{1}) && isfield(diminfo{1}, 'models')
    for i = 1:numel(diminfo{1}.models), G = G + diminfo{1}.models{i}.count; end
end
M = reshape(v, [rowLen nRows]);
if G > 0 && rowLen == G
    cdata = M;                 % [grayordinates x maps]
elseif G > 0 && nRows == G
    cdata = M';                % transpose to [grayordinates x maps]
else
    % Dense dimension not size-identified (e.g. dconn): keep file order,
    % assume dimension-1 (along rows) is the dense dimension.
    cdata = M;
end

cii.cdata   = cdata;
cii.intent  = intent;
cii.diminfo = diminfo;
cii.hdr     = hdr;

end % main


% =========================================================================
function prec = local_nifti_precision(dt)
switch dt
    case 2,    prec = 'uint8';
    case 4,    prec = 'int16';
    case 8,    prec = 'int32';
    case 16,   prec = 'single';
    case 64,   prec = 'double';
    case 256,  prec = 'int8';
    case 512,  prec = 'uint16';
    case 768,  prec = 'uint32';
    case 1024, prec = 'int64';
    case 1280, prec = 'uint64';
    otherwise
        error('canlab_read_cifti:datatype', 'Unsupported NIfTI datatype %d', dt);
end
end


% =========================================================================
function [diminfo, intent] = local_parse_cifti_xml(xml, intent_code, intent_name)
% Returns diminfo{1} (dense) and diminfo{2} (maps/series/labels).

switch intent_code
    case 3000, intent = 'dconn';
    case 3002, intent = 'dtseries';
    case 3006, intent = 'dscalar';
    case 3007, intent = 'dlabel';
    otherwise
        if     contains(intent_name, 'Scalar'), intent = 'dscalar';
        elseif contains(intent_name, 'Series'), intent = 'dtseries';
        elseif contains(intent_name, 'Label'),  intent = 'dlabel';
        elseif contains(intent_name, 'Dense'),  intent = 'dconn';
        else,  intent = 'unknown';
        end
end

% Match either a self-closing <MatrixIndicesMap .../> (e.g. SERIES maps carry
% all their info in attributes) or a paired <MatrixIndicesMap ...>...</...> block.
mims = regexp(xml, '(?s)<MatrixIndicesMap\s[^>]*?/>|<MatrixIndicesMap\s[^>]*?>.*?</MatrixIndicesMap>', 'match');

infos   = {};   % parsed info structs
applies = {};   % AppliesToMatrixDimension vectors
isdense = [];   % logical: this map is a BRAIN_MODELS (dense) dimension

for k = 1:numel(mims)
    mim = mims{k};
    appliesTo = local_attr(mim, 'AppliesToMatrixDimension'); % "0" / "1" / "0,1"
    maptype   = local_attr(mim, 'IndicesMapToDataType');

    info = struct();
    dense = false;
    if contains(maptype, 'BRAIN_MODELS')
        info = local_parse_brain_models(mim);
        dense = true;
    elseif contains(maptype, 'SCALARS')
        info.type = 'scalars';
        info.maps = local_parse_named_maps(mim, false);
        info.length = numel(info.maps);
    elseif contains(maptype, 'LABELS')
        info.type = 'labels';
        info.maps = local_parse_named_maps(mim, true);
        info.length = numel(info.maps);
    elseif contains(maptype, 'SERIES')
        info.type = 'series';
        info.seriesStart    = str2double(local_attr(mim, 'SeriesStart'));
        info.seriesStep     = str2double(local_attr(mim, 'SeriesStep'));
        info.seriesUnit     = local_attr(mim, 'SeriesUnit');
        info.seriesExponent = str2double(local_attr(mim, 'SeriesExponent'));
        info.length = str2double(local_attr(mim, 'NumberOfSeriesPoints'));
    elseif contains(maptype, 'PARCELS')
        info.type = 'parcels';   % not fully expanded in v1
        info.raw  = mim;
    else
        info.type = 'unknown';
    end

    infos{end+1}   = info;            %#ok<AGROW>
    applies{end+1} = sscanf(appliesTo, '%d,')'; %#ok<AGROW>
    isdense(end+1) = dense;           %#ok<AGROW>
end

% Order so diminfo{1} is the dense (grayordinate) dimension and diminfo{2} the
% maps/series/labels dimension -- matches cifti-matlab and get_cifti_data, which
% index diminfo{1}.models. (CIFTI AppliesToMatrixDimension numbering is reversed
% relative to MATLAB row/col, so we key off content, not the attribute.)
diminfo = {[], []};
denseIdx = find(isdense, 1);
mapsIdx  = find(~isdense, 1);
if ~isempty(denseIdx), diminfo{1} = infos{denseIdx}; end
if ~isempty(mapsIdx),  diminfo{2} = infos{mapsIdx};  end
if isempty(denseIdx)            % e.g. dconn: both dimensions are dense
    [~, order] = sort(cellfun(@(a) a(1), applies));
    diminfo = infos(order);
end
end


% =========================================================================
function info = local_parse_brain_models(mim)
info = struct('type', 'dense', 'length', 0, 'models', {{}}, 'vol', []);

% Volume block (affine + dims), if any
vol = regexp(mim, '(?s)<Volume\s.*?</Volume>', 'match', 'once');
if ~isempty(vol)
    vd = local_attr(vol, 'VolumeDimensions');
    dims = sscanf(vd, '%d,')';
    Telem = regexp(vol, '(?s)<TransformationMatrixVoxelIndicesIJKtoXYZ[^>]*>.*?</TransformationMatrixVoxelIndicesIJKtoXYZ>', ...
        'match', 'once');
    T = regexp(Telem, '(?s)>(.*?)<', 'tokens', 'once');
    sform = eye(4);
    if ~isempty(T)
        mexp = str2double(local_attr(Telem, 'MeterExponent'));
        if isnan(mexp), mexp = -3; end
        nums = sscanf(T{1}, '%f');
        if numel(nums) >= 16
            sform = reshape(nums(1:16), 4, 4)';        % row-major -> 4x4
            sform(1:3, :) = sform(1:3, :) * 10^(mexp + 3);  % normalize to mm
        end
    end
    info.vol = struct('dims', dims, 'sform', sform);
end

bms = regexp(mim, '(?s)<BrainModel\s.*?</BrainModel>', 'match');
total = 0;
for i = 1:numel(bms)
    bm = bms{i};
    m = struct();
    m.struct = regexprep(local_attr(bm, 'BrainStructure'), '^CIFTI_STRUCTURE_', '');
    mt = local_attr(bm, 'ModelType');
    if contains(mt, 'SURFACE'), m.type = 'surf'; else, m.type = 'vox'; end
    m.start   = str2double(local_attr(bm, 'IndexOffset')) + 1;  % 1-based
    m.count   = str2double(local_attr(bm, 'IndexCount'));
    m.numvert = str2double(local_attr(bm, 'SurfaceNumberOfVertices'));
    m.vertlist = [];
    m.voxlist  = [];
    if strcmp(m.type, 'surf')
        vi = regexp(bm, '(?s)<VertexIndices[^>]*>(.*?)</VertexIndices>', 'tokens', 'once');
        if ~isempty(vi), m.vertlist = sscanf(vi{1}, '%d')'; end   % 0-based
    else
        vx = regexp(bm, '(?s)<VoxelIndicesIJK[^>]*>(.*?)</VoxelIndicesIJK>', 'tokens', 'once');
        if ~isempty(vx)
            ijk = sscanf(vx{1}, '%d');
            m.voxlist = reshape(ijk, 3, [])';   % Nx3 IJK (0-based)
            m.voxlist = m.voxlist';             % 3xN to match cifti-matlab
        end
    end
    info.models{end+1} = m;
    total = total + m.count;
end
info.length = total;
end


% =========================================================================
function maps = local_parse_named_maps(mim, isLabel)
maps = struct('name', {}, 'table', {});
nm = regexp(mim, '(?s)<NamedMap>(.*?)</NamedMap>', 'match');
for i = 1:numel(nm)
    blk = nm{i};
    nameTok = regexp(blk, '(?s)<MapName>(.*?)</MapName>', 'tokens', 'once');
    name = '';
    if ~isempty(nameTok)
        name = strtrim(nameTok{1});
        name = regexprep(name, '<!\[CDATA\[(.*?)\]\]>', '$1');
    end
    tbl = [];
    if isLabel
        tbl = local_parse_cifti_labeltable(blk);
    end
    maps(end+1) = struct('name', name, 'table', tbl); %#ok<AGROW>
end
end


% =========================================================================
function tbl = local_parse_cifti_labeltable(blk)
tbl = struct('key', {}, 'name', {}, 'rgba', {});
lt = regexp(blk, '(?s)<LabelTable>(.*?)</LabelTable>', 'tokens', 'once');
if isempty(lt), return; end
items = regexp(lt{1}, '(?s)<Label\s.*?</Label>', 'match');
for i = 1:numel(items)
    it = items{i};
    key = str2double(local_attr(it, 'Key'));
    r = str2double(local_attr(it, 'Red'));
    g = str2double(local_attr(it, 'Green'));
    b = str2double(local_attr(it, 'Blue'));
    a = str2double(local_attr(it, 'Alpha'));
    nmt = regexp(it, '(?s)<Label\s[^>]*>(.*?)</Label>', 'tokens', 'once');
    name = '';
    if ~isempty(nmt)
        name = strtrim(nmt{1});
        name = regexprep(name, '<!\[CDATA\[(.*?)\]\]>', '$1');
    end
    tbl(end+1) = struct('key', key, 'name', name, 'rgba', [r g b a]); %#ok<AGROW>
end
end


% =========================================================================
function val = local_attr(str, name)
tok = regexp(str, [name '\s*=\s*"([^"]*)"'], 'tokens', 'once');
if isempty(tok), val = ''; else, val = tok{1}; end
end
