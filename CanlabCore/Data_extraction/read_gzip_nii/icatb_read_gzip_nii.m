function varargout = icatb_read_gzip_nii(fileIn, varargin)
%% Read data from gzip nifti file (*.nii.gz) directly using GZIPInputstream.
%
% Inputs:
%
% 1. fileIn - Nifti gzip file name
% 2. varargin - Variable arguments passed in pairs
%
%   a. read_hdr_only - Optional argument. By default, both header and data information is returned. If set to 1, only necessary fields in header
%   are returned (slices and timeNo variables are ignored).
%   b. slices - Axial slice/slices. You could enter row vector of slices. By
%   default, all slices are loaded.
%   c. timepoints - Timepoints of interest. By default, all timepoints are
%   returned.
%   d.mask - Mask of interest. Specify full file path or use boolean mask
%   or indices.
%   e. use_spm -  Gzip files are un-zipped and spm volume functions are used to read nifti
%   files
%   f. buffer_size - Default buffer size is set to 2^14 bytes (memory efficient) when reading
%   using GZIPInputstream. Max buffer size you could specify is 2^31-1.
%
% Outputs:
% [hdr, V] - if read_hdr_only == 1. hdr is header and V is spm volume
% structure.
% [data, hdr, V] if read_hdr_only == 0. Data is 4D array (x,y,z,t) or 2D
% array (x*y*z, t) if mask is specified, hdr is header and V os spm volume structure.
%
%


NIFTI_GZ = 0;

GZIPINFO.isLargeFile = 1;
GZIPINFO.buffer_size = 2^14;


%% Parse inputs
read_hdr_only = 0;
buffer_size = 2^14;
isLargeFile = 0;
try
    buffer_size = GZIPINFO.buffer_size;
catch
end

try
    isLargeFile = GZIPINFO.isLargeFile;
catch
end

isSPMRead = NIFTI_GZ;
if (isempty(isSPMRead))
    isSPMRead = 0;
end
slices = [];
timeNo = [];

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'read_hdr_only'))
        read_hdr_only = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'slices'))
        slices = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'timepoints') || strcmpi(varargin{n}, 'time'))
        timeNo = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'mask'))
        mask = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'use_spm'))
        isSPMRead = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'buffer_size'))
        buffer_size = varargin{n + 1};
    end
end

if (nargout > 3)
    error('Max number of output arguments can be returned is 3');
end

if (read_hdr_only)
    if (nargout > 2)
        error('Max number of output arguments can be returned is 2 with read header only');
    end
end

if (~usejava('jvm'))
    error('Requires jvm to read gzip nifti (*.nii.gz) file');
end

jar_p = fileparts(which('icatb_read_gzip_nii.m'));

if (~isdeployed)
    javaaddpath(fullfile(jar_p, 'icatb_gz.jar'));
end

import icatb_gz.*;

if (~isSPMRead)
    try
        zclass = icatb_gz.read_gzip();
    catch
        isSPMRead = 1;
    end
end

nBytesToRead = 348;
doSwapBytes = 0;

%% Use GZIPInputstream
if (read_hdr_only)
    fid  = java.io.FileInputStream(fileIn);
    zid  = java.util.zip.GZIPInputStream(fid);
    if (nargout == 1)
        zid.skip(40);
        if (~isSPMRead)
            inByteArray = zclass.read(zid, 16);
        else
            inByteArray = readInitialBytes(zid, 16);
        end
    else
        if (~isSPMRead)
            inByteArray = zclass.read(zid, nBytesToRead);
        else
            inByteArray = readInitialBytes(zid, nBytesToRead);
        end
    end
    try
        zid.close();
        fid.close();
    catch
    end
    hdr = read_hdr(inByteArray);
    if (nargout > 1)
        gzV = getVol(fileIn, hdr);
    end
    
else
    
    if (~isSPMRead)
        [inByteArray, hdr, doSwapBytes] = readGzip(fileIn, buffer_size, isLargeFile);
        gzV = getVol(fileIn, hdr);
    else
        % Uncompress file and read header
        [hdr, orig_gzfn, gzV, XYZ] = getSpmVol(fileIn, timeNo);
    end
    
end

if (read_hdr_only)
    
    varargout{1} = hdr;
    if (nargout == 2)
        varargout{2} = gzV;
    end
    
    return;
    
end

dims = hdr.dim(2:5);

if (~exist('slices', 'var') || isempty(slices))
    slices = 1:dims(3);
end

if (~exist('timeNo', 'var') || isempty(timeNo))
    timeNo = 1:dims(4);
end


byteA = double(hdr.bitpix / 8);

convertDataTo2D = 1;

%% Mask
if (exist('mask', 'var') && ~isempty(mask))
    if (ischar(mask))
        mask = spm_read_vols(spm_vol(mask));
        mask(isfinite(mask)==0) = 0;
    end
    if (numel(mask) ~= prod(dims(1:3)))
        maskF = zeros(dims(1:3));
        maskF(mask) = 1;
    else
        maskF = mask;
    end
    maskF = reshape(maskF, dims(1:3));
else
    convertDataTo2D = 0;
    maskF = (ones(dims(1:3)) ~= 0);
end

mask_inds = find(maskF(:, :, slices) ~= 0);

if (~convertDataTo2D)
    data = zeros(dims(1), dims(2), length(slices), length(timeNo));
else
    data = zeros(length(mask_inds), length(timeNo));
end


for t = 1:length(timeNo)
    endT = 0;
    for z = 1:length(slices)
        startT = endT  + 1;
        
        if (~isSPMRead)
            
            if (convertDataTo2D)
                mInds = find(maskF(:, :, slices(z)) ~= 0);
            end
            
            position =  double(byteA*((timeNo(t)-1)*prod(dims(1:3)) + (slices(z)-1)*prod(dims(1:2)))) + double(hdr.vox_offset);
            
            tmp = inByteArray(position + 1: position + dims(1)*dims(2)*byteA);
            
            if (~convertDataTo2D)
                if (~doSwapBytes)
                    tmp = double(typecast(tmp, hdr.precision));
                else
                    tmp = double(swapbytes(typecast(tmp, hdr.precision)));
                end
                tmpDat = reshape(tmp, [dims(1), dims(2)]);
            else
                if (~doSwapBytes)
                    tmpDat = double(typecast(tmp, hdr.precision));
                else
                    tmpDat = double(swapbytes(typecast(tmp, hdr.precision)));
                end
                tmpDat = tmpDat(mInds);
            end
            
            endT = endT + length(tmpDat);
            
        else
            
            mInds = find(maskF(:, :, slices(z)) ~= 0);
            tmpa = squeeze(XYZ(1, :, :, slices(z)));
            tmpb = squeeze(XYZ(2, :, :, slices(z)));
            tmpa = tmpa(mInds);
            tmpb = tmpb(mInds);
            tmpDat = [];
            if (~isempty(mInds))
                tmpDat = spm_sample_vol(gzV(t), tmpa, tmpb, slices(z)*ones(size(mInds)), 0);
            end
            endT = endT + length(tmpDat);
            if (~convertDataTo2D)
                tmpDat = reshape(tmpDat, dims(1), dims(2));
            end
            
        end
        
        clear tmp;
        
        if (~convertDataTo2D)
            data(:, :, z, t) = tmpDat;
        else
            data(startT:endT, t) = tmpDat;
        end
        
    end
    
end


if (~isSPMRead)
    if (hdr.scl_slope ~=0)
        data = data*hdr.scl_slope + hdr.scl_inter;
    end
end

data (isfinite(data) == 0) = 0;

varargout{1} = data;
varargout{2} = hdr;

if (nargout == 3)
    varargout{3} = gzV(1);
end

if (exist('orig_gzfn', 'var'))
    doCleanUpFiles(orig_gzfn);
end



function [hdr, doSwapBytes] = read_hdr(inByteArray)
%% Read header

doSwapBytes = 0;
corr_types = {'uint8', 'int16', 'int32', 'single', 'double', 'int8', 'uint16', 'uint32'};
codes   = [2, 4, 8, 16, 64, 256, 512, 768];

% necessary fields
if (length(inByteArray) == 16)
    dim = typecast(inByteArray, 'int16');
    if (dim(1) < 1 || dim(1) > 7)
        doSwapBytes = 1;
        dim = swapbytes(typecast(inByteArray, 'int16'));
    end
    hdr.dim = dim(:)';
    return;
end


chk = typecast(inByteArray(41:56), 'int16');
if (chk(1) < 1 || chk(1) > 7)
    doSwapBytes = 1;
end

fieldsIn  = {   'dim',    'int16',  (41:56);
    'datatype',           'int16',  (71:72);
    'bitpix',             'int16',  (73:74);
    'vox_offset',         'single', (109:112);
    'pixdim',             'single', (77:108);
    'scl_slope',          'single', (113:116);
    'scl_inter',          'single', (117:120);
    'qform_code',         'int16',  (253:254);
    'sform_code',         'int16',  (255:256);
    'quatern_b',          'single', (257:260);
    'quatern_c',          'single', (261:264);
    'quatern_d',          'single', (265:268);
    'qoffset_x',          'single', (269:272);
    'qoffset_y',          'single', (273:276);
    'qoffset_z',          'single', (277:280);
    'sizeof_hdr',         'int32',  (1:4);
    'data_type',          'char',   (5:14);
    'db_name',            'char',   (15:32);
    'extents',            'int32',  (33:36);
    'session_error',      'int16',  (37:38);
    'regular',            'char',   (39:39);
    'dim_info',           'uint8',   (40:40);
    'intent_p1',          'single', (57:60);
    'intent_p2',          'single', (61:64);
    'intent_p3',          'single', (65:68);
    'intent_code',        'int16',  (69:70);
    'slice_start',        'int16',  (75:76);
    'slice_end',          'int16',  (121:122);
    'slice_code',         'uint8',   (123:123);
    'xyzt_units',         'uint8',   (124:124);
    'cal_max',            'single', (125:128);
    'cal_min',            'single', (129:132);
    'slice_duration',     'single', (133:136);
    'toffset',            'single', (137:140);
    'glmax',              'single', (141:144);
    'glmin',              'single', (145:148);
    'descrip',            'char',   (149:228);
    'aux_file',           'char',   (229:252);
    'srow_x',             'single', (281:296);
    'srow_y',             'single', (297:312);
    'srow_z',             'single', (313:328);
    'intent_name',        'char',   (329:344);
    'magic',              'char',   (345:348)};

for nF = 1:size(fieldsIn, 1)
    if (strcmpi(fieldsIn{nF, 2}, 'char'))
        tmp = char(inByteArray(fieldsIn{nF, 3}));
        if (size(tmp, 1) > 1)
            tmp = tmp';
        end
    else
        if (~doSwapBytes)
            tmp = typecast(inByteArray(fieldsIn{nF, 3}), fieldsIn{nF, 2});
        else
            tmp = swapbytes(typecast(inByteArray(fieldsIn{nF, 3}), fieldsIn{nF, 2}));
        end
    end
    if (isnumeric(tmp))
        tmp = double(tmp(:)');
    end
    hdr.(fieldsIn{nF, 1}) = tmp;
end

hdr.dim = hdr.dim(:)';
hdr.pixdim = hdr.pixdim(:)';

hdr.precision = [];
chk = find(codes == hdr.datatype);
if (~isempty(chk))
    hdr.precision = corr_types{chk};
end

if isempty(hdr.precision)
    error('Unknown data type');
end


function [inByteArray, hdr, doSwapBytes] = readGzip(fileIn, buffer_size, isLargeFile)
%% Read large data in chunks
%

import icatb_gz.*;
zclass = icatb_gz.read_gzip();
hdrSize = 348;

tmpFileInfo = dir(fileIn);
tmpGBytes = tmpFileInfo(1).bytes/1024/1024/1024;

% if (tmpGBytes > 0.5)
%     isLargeFile = 1;
% end


if (isLargeFile)
    
    fid  = java.io.FileInputStream(fileIn);
    zid  = java.util.zip.GZIPInputStream(fid);
    
    try
        
        hdrBytes = zclass.read(zid, hdrSize);
        
        [hdr, doSwapBytes] = read_hdr(typecast(hdrBytes, 'uint8'));
        
        voxel_offset = hdr.vox_offset;
        bytesA = hdr.bitpix/8;
        byte_array_size = voxel_offset + prod(hdr.dim(2:5))*bytesA;
        inByteArray = zeros(byte_array_size, 1, 'int8');
        inByteArray(1:hdrSize) = hdrBytes;
        zid.close();
        fid.close();
        
        fid  = java.io.FileInputStream(fileIn);
        zid  = java.util.zip.GZIPInputStream(fid);
        zid.skip(voxel_offset);
        count = voxel_offset;
        while 1
            
            count = count + 1;
            
            tmp = zclass.read(zid, buffer_size);
            
            len = length(tmp);
            
            if (len == 0)
                break;
            end
            
            endCount = count + len - 1;
            inByteArray(count:endCount) = tmp;
            count = endCount;
            
        end
        
        inByteArray = typecast(inByteArray, 'uint8');
        
        zid.close();
        fid.close();
        
    catch
        
        
        zid.close();
        fid.close();
        
    end
    
    
else
    
    inByteArray = zclass.read(fileIn, buffer_size);
    inByteArray = typecast(inByteArray, 'uint8');
    [hdr, doSwapBytes] = read_hdr(inByteArray);
    
end

function V = getVol(fname, hdr)
%% Create spm volume structure from header

fname_i = fname;
[pathstr, fn, extn] = fileparts(fname);

extn = strrep(lower(extn), '.gz', '');
fname = fullfile(pathstr, [fn, extn]);

mat = decode_qform0(hdr);

V   = struct('fname', fname,...
    'dim',   hdr.dim(2:4),...
    'dt',   [hdr.datatype, 0],...
    'pinfo', [hdr.scl_slope hdr.scl_inter hdr.vox_offset]',...
    'mat',   mat,...
    'n',     1,...
    'descrip', 'NIFTI');

h = nifti;
h.mat = V.mat;
h.mat0 = V.mat;
fp  = fopen(fname_i, 'r', 'native');
[~,~,mach] = fopen(fp);
fclose(fp);
dat = file_array(V.fname, hdr.dim(2:5),[hdr.datatype, strfind(mach,'be')], hdr.vox_offset, hdr.scl_slope, hdr.scl_inter);
h.dat = dat;
h.diminfo.slice = 3;
h.diminfo.slice_time.code = hdr.slice_code;
h.diminfo.slice_time.start = hdr.slice_start+1;
h.diminfo.slice_time.end = hdr.slice_end+1;
h.timing.toffset = hdr.toffset;
h.timing.tspace = hdr.pixdim(5);
matcodes = {'UNKNOWN', 'Scanner', 'Aligned', 'Talairach', 'MNI152'};
h.mat_intent = matcodes{hdr.sform_code+1};
h.mat0_intent = matcodes{hdr.qform_code+1};
V.private = h;


function M = decode_qform0(hdr)
% Decode qform info from NIFTI-1 headers.
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: decode_qform0.m 3131 2009-05-18 15:54:10Z guillaume $


dim    = double(hdr.dim);
pixdim = double(hdr.pixdim);
if ~isfield(hdr,'magic') || hdr.qform_code <= 0,
    flp = spm_flip_analyze_images;
    %disp('------------------------------------------------------');
    %disp('The images are in a form whereby it is not possible to');
    %disp('tell the left and right sides of the brain apart.');
    %if flp,
    %    disp('They are assumed to be stored left-handed.');
    %else
    %    disp('They are assumed to be stored right-handed.');
    %end;
    %disp('------------------------------------------------------');
    
    %R     = eye(4);
    n      = min(dim(1),3);
    vox    = [pixdim(2:(n+1)) ones(1,3-n)];
    
    if ~isfield(hdr,'origin') || ~any(hdr.origin(1:3)),
        origin = (dim(2:4)+1)/2;
    else
        origin = double(hdr.origin(1:3));
    end;
    off     = -vox.*origin;
    M       = [vox(1) 0 0 off(1) ; 0 vox(2) 0 off(2) ; 0 0 vox(3) off(3) ; 0 0 0 1];
    
    % Stuff for default orientations
    if flp, M = diag([-1 1 1 1])*M; end;
else
    
    % Rotations from quaternions
    R = Q2M(double([hdr.quatern_b hdr.quatern_c hdr.quatern_d]));
    
    % Translations
    T = [eye(4,3) double([hdr.qoffset_x hdr.qoffset_y hdr.qoffset_z 1]')];
    
    % Zooms.  Note that flips are derived from the first
    % element of pixdim, which is normally unused.
    n = min(dim(1),3);
    Z = [pixdim(2:(n+1)) ones(1,4-n)];
    Z(Z<0) = 1;
    if pixdim(1)<0, Z(3) = -Z(3); end;
    Z = diag(Z);
    
    M = T*R*Z;
    
    % Convert from first voxel at [1,1,1]
    % to first voxel at [0,0,0]
    M = M * [eye(4,3) [-1 -1 -1 1]'];
end;
return;



function M = Q2M(Q)
% Generate a rotation matrix from a quaternion xi+yj+zk+w,
% where Q = [x y z], and w = 1-x^2-y^2-z^2.
% See: http://skal.planet-d.net/demo/matrixfaq.htm
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: Q2M.m 1143 2008-02-07 19:33:33Z spm $


Q = Q(1:3); % Assume rigid body
w = sqrt(1 - sum(Q.^2));
x = Q(1); y = Q(2); z = Q(3);
if w<1e-7,
    w = 1/sqrt(x*x+y*y+z*z);
    x = x*w;
    y = y*w;
    z = z*w;
    w = 0;
end;
xx = x*x; yy = y*y; zz = z*z; ww = w*w;
xy = x*y; xz = x*z; xw = x*w;
yz = y*z; yw = y*w; zw = z*w;
M = [...
    (xx-yy-zz+ww)      2*(xy-zw)      2*(xz+yw) 0
    2*(xy+zw) (-xx+yy-zz+ww)      2*(yz-xw) 0
    2*(xz-yw)      2*(yz+xw) (-xx-yy+zz+ww) 0
    0              0              0  1];
return;


function varargout = getSpmVol(fileIn, timeNo)
%% Unzip file and read info
%

gzfn = gunzip (fileIn, tempdir);
orig_gzfn = char(gzfn);

if (nargout > 2)
    gzfn = rename_4d_file(orig_gzfn);
    if (~exist('timeNo', 'var') || isempty(timeNo))
        timeNo = (1:size(gzfn, 1));
    end
    gzV = spm_vol(deblank(gzfn(timeNo, :)));
    hdr = gzV(1).private.hdr;
    
    hdr.dim = double(hdr.dim);
    
else
    
    hdr = icatb_read_hdr(deblank(orig_gzfn));
    hdr.dim = double(hdr.dime.dim);
    
end

varargout{1} = hdr;
varargout{2} = orig_gzfn;

if (nargout >= 3)
    varargout{3} = gzV;
end

if (nargout == 4)
    XYZ = icatb_get_voxel_coords(gzV(1).dim(1:3));
    XYZ = reshape(XYZ, [3, gzV(1).dim(1:3)]);
    varargout{4} = XYZ;
end


function doCleanUpFiles(orig_gzfn)

try
    for nF = 1:size(orig_gzfn, 1)
        delete(deblank(orig_gzfn(nF, :)));
    end
catch
end

function inByteArray = readInitialBytes(zid, nBytesToRead)
%% Use read method in matlab
%

if (isfinite(nBytesToRead))
    inByteArray = zeros(nBytesToRead, 1, 'uint8');
else
    inByteArray = [];
    inByteArray = uint8(inByteArray);
end

count = 0;
tmp = 0;
while ((tmp ~= -1) && (count < nBytesToRead))
    count = count + 1;
    tmp = uint8(read(zid));
    inByteArray(count) = tmp;
end





function [hdr,otherendian] = icatb_read_hdr(fname)
% Read (SPM customised) Analyze header
% fname       - .hdr filename
% hdr         - structure containing Analyze header
% otherendian - byte swapping necessary flag
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_read_hdr.m 112 2005-05-04 18:20:52Z john $

otherendian = 0;
fname = deblank(fname);
if (strcmpi(fname(end-2:end), '.gz'))
    hdr = icatb_read_gzip_nii(fname, 'read_hdr_only', 1);
    return;
end


fid         = fopen(fname,'r','native');
%otherendian = 0;
if (fid > 0)
    dime = read_dime(fid);
    if dime.dim(1)<0 | dime.dim(1)>15, % Appears to be other-endian
        % Re-open other-endian
        fclose(fid);
        if spm_platform('bigend'), fid = fopen(fname,'r','ieee-le');
        else,                      fid = fopen(fname,'r','ieee-be'); end;
        otherendian = 1;
        dime = read_dime(fid);
    end;
    %hk       = read_hk(fid);
    %hist     = read_hist(fid);
    %hdr.hk   = hk;
    hdr.dime = dime;
    %hdr.hist = hist;
    
    % SPM specific bit - unused
    %if hdr.hk.sizeof_hdr > 348,
    %	spmf = read_spmf(fid,dime.dim(5));
    %	if ~isempty(spmf),
    %		hdr.spmf = spmf;
    %	end;
    %end;
    
    fclose(fid);
else,
    hdr = [];
    otherendian = NaN;
    %error(['Problem opening header file (' fopen(fid) ').']);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function hk = read_hk(fid)
% read (struct) header_key
%-----------------------------------------------------------------------
fseek(fid,0,'bof');
hk.sizeof_hdr 		= fread(fid,1,'int32');
hk.data_type  		= mysetstr(fread(fid,10,'uchar'))';
hk.db_name    		= mysetstr(fread(fid,18,'uchar'))';
hk.extents    		= fread(fid,1,'int32');
hk.session_error	= fread(fid,1,'int16');
hk.regular			= mysetstr(fread(fid,1,'uchar'))';
hk.hkey_un0			= mysetstr(fread(fid,1,'uchar'))';
if isempty(hk.hkey_un0), error(['Problem reading "hk" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function dime = read_dime(fid)
% read (struct) image_dimension
%-----------------------------------------------------------------------
fseek(fid,40,'bof');
dime.dim		= fread(fid,8,'int16')';
%dime.vox_units	= mysetstr(fread(fid,4,'uchar'))';
%dime.cal_units	= mysetstr(fread(fid,8,'uchar'))';
%dime.unused1	= fread(fid,1,'int16');
%dime.datatype	= fread(fid,1,'int16');
%dime.bitpix		= fread(fid,1,'int16');
%dime.dim_un0	= fread(fid,1,'int16');
%dime.pixdim		= fread(fid,8,'float')';
%dime.vox_offset	= fread(fid,1,'float');
%dime.funused1	= fread(fid,1,'float');
%dime.funused2	= fread(fid,1,'float');
%dime.funused3	= fread(fid,1,'float');
%dime.cal_max	= fread(fid,1,'float');
%dime.cal_min	= fread(fid,1,'float');
%dime.compressed	= fread(fid,1,'int32');
%dime.verified	= fread(fid,1,'int32');
dime.glmax		= fread(fid,1,'int32');
dime.glmin		= fread(fid,1,'int32');
if isempty(dime.glmin), error(['Problem reading "dime" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function hist = read_hist(fid)
% read (struct) data_history
%-----------------------------------------------------------------------
fseek(fid,148,'bof');
hist.descrip	= mysetstr(fread(fid,80,'uchar'))';
hist.aux_file	= mysetstr(fread(fid,24,'uchar'))';
hist.orient		= fread(fid,1,'uchar');
hist.origin		= fread(fid,5,'int16')';
hist.generated	= mysetstr(fread(fid,10,'uchar'))';
hist.scannum	= mysetstr(fread(fid,10,'uchar'))';
hist.patient_id	= mysetstr(fread(fid,10,'uchar'))';
hist.exp_date	= mysetstr(fread(fid,10,'uchar'))';
hist.exp_time	= mysetstr(fread(fid,10,'uchar'))';
hist.hist_un0	= mysetstr(fread(fid,3,'uchar'))';
hist.views		= fread(fid,1,'int32');
hist.vols_added	= fread(fid,1,'int32');
hist.start_field= fread(fid,1,'int32');
hist.field_skip	= fread(fid,1,'int32');
hist.omax		= fread(fid,1,'int32');
hist.omin		= fread(fid,1,'int32');
hist.smax		= fread(fid,1,'int32');
hist.smin		= fread(fid,1,'int32');
if isempty(hist.smin), error(['Problem reading "hist" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function spmf = read_spmf(fid,n)
% Read SPM specific fields
% This bit may be used in the future for extending the Analyze header.

fseek(fid,348,'bof');
mgc = fread(fid,1,'int32');    % Magic number
if mgc ~= 20020417, spmf = []; return; end;

for j=1:n,
    spmf(j).mat    = fread(fid,16,'double'); % Orientation information
    spmf(j).unused = fread(fid,384,'uchar'); % Extra unused stuff
    if length(spmf(j).unused)<384,
        error(['Problem reading "spmf" of header file (' fopen(fid) ').']);
    end;
    spmf(j).mat = reshape(spmf(j).mat,[4 4]);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function out = mysetstr(in)
tmp = find(in == 0);
tmp = min([min(tmp) length(in)]);
out = setstr([in(1:tmp)' zeros(1,length(in)-(tmp))])';
return;
%_______________________________________________________________________
%_______________________________________________________________________




function files = rename_4d_file(files, fileNumber)
%For nifti and analyze files rename the files by adding a number at the
%end

if ~exist('fileNumber', 'var')
    fileNumber = [];
end

if ~isempty(files)
    
    % Do pattern match to check IMG or NII files
    files = cellstr(files);
    
    % MATCH IMG OR NII
    checkNII = regexpi(files, '\.nii$|\.gz$');
    checkIMG = regexpi(files, '\.img$');
    
    good_inds1 = icatb_good_cells(checkNII);
    good_inds2 = icatb_good_cells(checkIMG);
    
    % Get good cells
    good_inds = (good_inds1 | good_inds2);
    good_inds = find(good_inds ~= 0);
    
    clear good_inds1 good_inds2 checkNII checkIMG;
    
    % If IMG or NII files exist check the headers of the image files
    if ~isempty(good_inds)
        imFiles = files(good_inds);
        if (~isempty(fileNumber)) && (length(fileNumber) == 1) && (fileNumber == 1)
            imFiles = strcat(imFiles, ',1');
            files(good_inds) = imFiles;
        else
            % Loop over img or nii files
            for nIm = 1:length(imFiles)
                currentFile = imFiles{nIm};
                
                try
                    numFiles = getDims(currentFile);
                catch
                    numFiles = 1;
                end
                
                if (isempty(fileNumber))
                    tempFileNum = (1:numFiles);
                else
                    tempFileNum = fileNumber;
                end
                
                tempFileNum(tempFileNum > numFiles) = [];
                
                if ~isempty(tempFileNum)
                    tempFiles = [repmat(currentFile, length(tempFileNum), 1), repmat(',', length(tempFileNum), 1), numberToString(tempFileNum(:))];
                else
                    tempFiles = '';
                end
                files{good_inds(nIm)} = tempFiles;
            end
            % End loop over img or nii files
            files = cellstr(char(files));
            ind = icatb_good_cells(files);
            files = files(ind);
        end
    end
    % End for checking IMG or NII files
    
    if ~isempty(files)
        files = char(files);
    end
    
end


function str = numberToString(nums)
% number to string. pad space after the number

try
    str = arrayfun(@num2str, nums, 'uniformoutput', false);
catch
    str = cell(length(nums), 1);
    for n = 1:length(nums)
        str{n} = num2str(nums(n));
    end
end

str = char(str);


function tp = getDims(currentFile)

if (~strcmpi(currentFile(end-2:end),'.gz'))
    ni = nifti(currentFile);
    tp = ni.dat.dim(4);
else
    ni = icatb_read_hdr(currentFile);
    tp = ni.dim(5);
end


function ind = icatb_good_cells(mycell)
%% Find good cells

if ~iscell(mycell)
    mycell = {mycell};
end

ind = cellfun('isempty', mycell);

% Good cells
ind = (ind == 0);

% Convert ind to row vector
ind = ind(:)';


function xyz = icatb_get_voxel_coords(dims)
%% Get voxel coords
%

xdim = dims(1);
ydim = dims(2);
zdim = dims(3);

[xords, yords] = ndgrid(1:xdim, 1:ydim);
xords = xords(:)'; yords = yords(:)';
CrPl    = 1:zdim;
zords   = CrPl(:)*ones(1,xdim*ydim);
xyz   = [repmat(xords,1,numel(CrPl)); ...
    repmat(yords,1,numel(CrPl)); ...
    reshape(zords',1,[])];