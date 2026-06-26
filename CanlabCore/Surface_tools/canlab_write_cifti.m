function canlab_write_cifti(filename, cii)
% Write a CIFTI-2 file (.dscalar/.dtseries/.dlabel.nii) natively in MATLAB.
%
% :Usage:
% ::
%     canlab_write_cifti(filename, cii)
%
% Writes a NIfTI-2 container with the CIFTI XML in an ecode==32 header
% extension, followed by the data matrix. The input struct matches
% canlab_read_cifti's output, so read->write->read round-trips exactly.
%
% Two modes:
%   * Round-trip: if `cii` carries the original `.xml` and `.hdr` (as returned
%     by canlab_read_cifti), the original XML and matrix layout are preserved
%     verbatim for a faithful re-write.
%   * Regenerate: otherwise the CIFTI XML is rebuilt from `cii.diminfo` (dense
%     BrainModels + a Scalars/Series/Labels map) in the standard layout.
%
% Depends only on core MATLAB -- no cifti-matlab, FieldTrip, or wb_command.
%
% :Inputs:
%
%   **filename:**
%        Char/string. Output path (e.g. '/tmp/out.dscalar.nii').
%
%   **cii:**
%        Struct (as from canlab_read_cifti) with at least:
%          .cdata    [nGrayordinates x nMaps]
%          .intent   'dscalar'|'dtseries'|'dlabel'|'dconn'
%          .diminfo  {dense, maps} as described in canlab_read_cifti
%        Optional .xml and .hdr trigger faithful round-trip mode.
%
% :Examples:
% ::
%     c = canlab_read_cifti(which('transcriptomic_gradients.dscalar.nii'));
%     c.cdata = zscore(c.cdata);                 % modify
%     canlab_write_cifti('/tmp/z.dscalar.nii', c);
%
% :See also: canlab_read_cifti, canlab_read_gifti, canlab_write_gifti

% ..
%    Part of the CANlab fmri_surface_data surface-data toolset. See
%    docs/fmri_surface_data_design_plan.md. No external toolbox required.
% ..

if nargin < 2 || ~isstruct(cii) || ~isfield(cii, 'cdata')
    error('canlab_write_cifti:input', 'Provide (filename, cii) with cii.cdata. See help.');
end
filename = char(filename);
cdata = cii.cdata;
G = size(cdata, 1);          % grayordinates
nMaps = size(cdata, 2);

% ---- XML: re-emit original or regenerate ----
reuse = isfield(cii, 'xml') && ~isempty(cii.xml) && isfield(cii, 'hdr') ...
        && isfield(cii.hdr, 'dim');
if reuse
    xml = cii.xml;
    rowLen = cii.hdr.dim(6);
    nRows  = cii.hdr.dim(7);
    % Reconstruct the on-disk vector in the file's original orientation
    if rowLen == G
        M = cdata;                 % cdata was stored as [rowLen x nRows]
    else
        M = cdata.';               % cdata was transposed on read
    end
    v = M(:);
    ndim = cii.hdr.dim(1);
else
    xml = local_build_cifti_xml(cii);
    rowLen = nMaps;                 % standard layout: maps fast, grayords slow
    nRows  = G;
    v = reshape(cdata.', [], 1);    % maps fastest within each grayordinate
    ndim = 6;
end

[intent_code, intent_name] = local_intent_codes(cii.intent);

% ---- Build header extension (ecode 32, padded so esize is a multiple of 16) ----
xmlbytes = uint8(xml);
raw_len = numel(xmlbytes);
esize = ceil((raw_len + 8) / 16) * 16;       % includes the 8-byte esize+ecode
pad = esize - 8 - raw_len;
vox_offset = 544 + esize;                     % 540 hdr + 4 extender + esize

% ---- NIfTI-2 dim vector ----
dim = ones(1, 8);
dim(1) = ndim;
dim(6) = rowLen;     % NIfTI dim[5]
dim(7) = nRows;      % NIfTI dim[6]

% ---- Write ----
fid = fopen(filename, 'w', 'l');
if fid < 0, error('canlab_write_cifti:open', 'Cannot open for writing: %s', filename); end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

hdrbuf = local_nifti2_header(dim, 16, 32, vox_offset, intent_code, intent_name);
fwrite(fid, hdrbuf, 'uint8');

fwrite(fid, uint8([1 0 0 0]), 'uint8');       % extension flag
fwrite(fid, esize, 'int32');
fwrite(fid, 32, 'int32');                     % ecode 32 = CIFTI
fwrite(fid, xmlbytes, 'uint8');
if pad > 0, fwrite(fid, zeros(1, pad, 'uint8'), 'uint8'); end

% Data, single precision (scl_slope written as 0 => values are literal)
here = ftell(fid);
if here < vox_offset, fwrite(fid, zeros(1, vox_offset - here, 'uint8'), 'uint8'); end
fwrite(fid, single(v), 'single');

end % main


% =========================================================================
function buf = local_nifti2_header(dim, datatype, bitpix, vox_offset, intent_code, intent_name)
buf = zeros(1, 540, 'uint8');
put = @(off, val, type) local_put(buf, off, val, type);

buf = put(0,  540, 'int32');
% magic "n+2\0" 0x0D 0x0A 0x1A 0x0A
buf(5:12) = uint8([110 43 50 0 13 10 26 10]);
buf = local_put(buf, 12, datatype, 'int16');
buf = local_put(buf, 14, bitpix,   'int16');
buf = local_put(buf, 16, int64(dim(:)'), 'int64');     % dim[8] @16
buf = local_put(buf, 104, ones(1, 8), 'double');       % pixdim[8]=1
buf = local_put(buf, 168, int64(vox_offset), 'int64'); % vox_offset
buf = local_put(buf, 176, 0, 'double');                % scl_slope=0 (no scaling)
buf = local_put(buf, 184, 0, 'double');                % scl_inter
buf = local_put(buf, 344, 0, 'int32');                 % qform_code
buf = local_put(buf, 348, 0, 'int32');                 % sform_code
buf = local_put(buf, 504, int32(intent_code), 'int32');
nm = uint8(intent_name); nm = nm(1:min(15, numel(nm)));
buf(509:509+numel(nm)-1) = nm;                          % intent_name @508 (16)
end


% =========================================================================
function buf = local_put(buf, off, val, type)
bytes = typecast(cast(val(:)', type), 'uint8');
buf(off+1:off+numel(bytes)) = bytes;
end


% =========================================================================
function [code, name] = local_intent_codes(intent)
switch intent
    case 'dscalar',  code = 3006; name = 'ConnDenseScalar';
    case 'dtseries', code = 3002; name = 'ConnDenseSeries';
    case 'dlabel',   code = 3007; name = 'ConnDenseLabel';
    case 'dconn',    code = 3000; name = 'ConnDense';
    otherwise,       code = 3006; name = 'ConnDenseScalar';
end
end


% =========================================================================
function xml = local_build_cifti_xml(cii)
% Regenerate CIFTI-2 XML from cii.diminfo (maps dim 0, dense BrainModels dim 1).
dense = cii.diminfo{1};
maps  = cii.diminfo{2};

s = sprintf('<?xml version="1.0" encoding="UTF-8"?>\n<CIFTI Version="2">\n<Matrix>\n');

% --- maps dimension (AppliesToMatrixDimension=0) ---
switch maps.type
    case 'scalars'
        s = [s sprintf('<MatrixIndicesMap AppliesToMatrixDimension="0" IndicesMapToDataType="CIFTI_INDEX_TYPE_SCALARS">\n')];
        for k = 1:numel(maps.maps)
            s = [s sprintf('<NamedMap><MapName>%s</MapName></NamedMap>\n', local_xesc(maps.maps(k).name))]; %#ok<AGROW>
        end
        s = [s sprintf('</MatrixIndicesMap>\n')];
    case 'labels'
        s = [s sprintf('<MatrixIndicesMap AppliesToMatrixDimension="0" IndicesMapToDataType="CIFTI_INDEX_TYPE_LABELS">\n')];
        for k = 1:numel(maps.maps)
            s = [s sprintf('<NamedMap><MapName>%s</MapName>\n<LabelTable>\n', local_xesc(maps.maps(k).name))]; %#ok<AGROW>
            t = maps.maps(k).table;
            for j = 1:numel(t)
                s = [s sprintf('<Label Key="%d" Red="%g" Green="%g" Blue="%g" Alpha="%g">%s</Label>\n', ...
                    t(j).key, t(j).rgba(1), t(j).rgba(2), t(j).rgba(3), t(j).rgba(4), local_xesc(t(j).name))]; %#ok<AGROW>
            end
            s = [s sprintf('</LabelTable>\n</NamedMap>\n')]; %#ok<AGROW>
        end
        s = [s sprintf('</MatrixIndicesMap>\n')];
    case 'series'
        s = [s sprintf(['<MatrixIndicesMap AppliesToMatrixDimension="0" IndicesMapToDataType="CIFTI_INDEX_TYPE_SERIES" ' ...
            'NumberOfSeriesPoints="%d" SeriesExponent="%d" SeriesStart="%g" SeriesStep="%g" SeriesUnit="%s"/>\n'], ...
            size(cii.cdata,2), local_def(maps,'seriesExponent',0), local_def(maps,'seriesStart',0), ...
            local_def(maps,'seriesStep',1), local_defstr(maps,'seriesUnit','SECOND'))];
    otherwise
        error('canlab_write_cifti:maptype', 'Cannot regenerate XML for maps type "%s". Re-write using the original .xml.', maps.type);
end

% --- dense dimension (AppliesToMatrixDimension=1) ---
s = [s sprintf('<MatrixIndicesMap AppliesToMatrixDimension="1" IndicesMapToDataType="CIFTI_INDEX_TYPE_BRAIN_MODELS">\n')];
if isfield(dense, 'vol') && ~isempty(dense.vol)
    A = dense.vol.sform;
    vals = A'; vals = vals(:)';   % row-major 16
    s = [s sprintf('<Volume VolumeDimensions="%d,%d,%d">\n', dense.vol.dims(1), dense.vol.dims(2), dense.vol.dims(3))];
    s = [s sprintf('<TransformationMatrixVoxelIndicesIJKtoXYZ MeterExponent="-3">')];
    s = [s sprintf('%.10g ', vals)];
    s = [s sprintf('</TransformationMatrixVoxelIndicesIJKtoXYZ>\n</Volume>\n')];
end
for i = 1:numel(dense.models)
    m = dense.models{i};
    if strcmp(m.type, 'surf')
        s = [s sprintf(['<BrainModel IndexOffset="%d" IndexCount="%d" ModelType="CIFTI_MODEL_TYPE_SURFACE" ' ...
            'BrainStructure="CIFTI_STRUCTURE_%s" SurfaceNumberOfVertices="%d">\n<VertexIndices>'], ...
            m.start-1, m.count, m.struct, m.numvert)]; %#ok<AGROW>
        s = [s sprintf('%d ', m.vertlist)];
        s = [s sprintf('</VertexIndices>\n</BrainModel>\n')]; %#ok<AGROW>
    else
        s = [s sprintf(['<BrainModel IndexOffset="%d" IndexCount="%d" ModelType="CIFTI_MODEL_TYPE_VOXELS" ' ...
            'BrainStructure="CIFTI_STRUCTURE_%s">\n<VoxelIndicesIJK>'], ...
            m.start-1, m.count, m.struct)]; %#ok<AGROW>
        ijk = m.voxlist;                  % 3xN
        s = [s sprintf('%d %d %d\n', ijk)];
        s = [s sprintf('</VoxelIndicesIJK>\n</BrainModel>\n')]; %#ok<AGROW>
    end
end
s = [s sprintf('</MatrixIndicesMap>\n</Matrix>\n</CIFTI>\n')];
xml = s;
end


function v = local_def(s, f, d), if isfield(s,f) && ~isnan(s.(f)), v = s.(f); else, v = d; end, end
function v = local_defstr(s, f, d), if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end, end


% =========================================================================
function out = local_xesc(str)
out = str;
out = strrep(out, '&', '&amp;');
out = strrep(out, '<', '&lt;');
out = strrep(out, '>', '&gt;');
end
