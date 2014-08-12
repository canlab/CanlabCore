function V = mask2struct(maskname,varargin)
% function V = mask2struct(maskname,crit_t,cl_size)
% inputs:
%	maskname: name of spmT, con, or filtered image
%	without .img extension, in single quotes
% optional inputs:
%	crit_t, cl_size: 	critical t and cluster size at which to mask
% 
% output: structure compatible with SPM viewing
% and with cluster definition algorithm tor_extract_rois
% to extract clusters:
% [clusters] = tor_extract_rois(maskname,V,V);
%
% to display:
% spm_image	(and choose anatomical)
% spm_orthviews('AddBlobs',1,V.XYZ,V.Z,V.mat)
% spm_orthviews('AddColouredBlobs',1,V.XYZ,V.Z,V.mat,[0 0 1])
%
% to overlay on Talairach atlas
% fixed_TSU(clusters)
% 
crit_t = 0;
cl_size = 0;
numClusters = NaN;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% input arguments
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if nargin > 1, crit_t = varargin{1};, end
if nargin > 2, cl_size = varargin{2};, end


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% read in the mask file
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
inV = spm_vol(maskname);
[vol,hdr] = readim2(maskname);
vol = spm_read_vols(inV);       % deal with SPM scaling factor!
vol = double(vol);


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% turn mask into 1s and 0s
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[maskedImage maskingImage] = maskImg(vol,crit_t,Inf);
	% if crit_t = 0
	%maskedImage is all positive values
	%maskingImage has all positive values = 1

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% make SPM.mat-like structure V
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	% Get XYZ pointlist from mask and cluster size mask if spec
	% ------------------------------------------------------
	if cl_size > 0
		[maskingImage,numClusters,XYZ] = clusterSizeMask(cl_size,maskingImage);
	else
		XYZ = mask2voxel(maskingImage);
		XYZ = XYZ';
	end

	% Get Z output - intensity values for sig. voxels
	% ------------------------------------------------------
	for i = 1:size(XYZ,2)
		% row is y, col is x
		Z(i) = vol(XYZ(1,i),XYZ(2,i),XYZ(3,i));
    end
    
    if isempty(XYZ)
        Z = [];
    end

	% get the header and volume information from the maskname
	% ------------------------------------------------------
	V = spm_vol(maskname);

	% adjust crit_t from 0 to 1 for compatibility with fixed_TSU (Talairach)
	% needs to be a non-zero value to display on Talairach
	% ------------------------------------------------------
	if crit_t == 0, crit_t = 1;, end

	% add other fields of V
	% ------------------------------------------------------
	V.maskingImage = maskingImage;
	V.cl_size = cl_size;
	V.crit_t = crit_t;
	V.numClusters = numClusters;
	V.XYZ = XYZ;
	V.Z = Z;
	V.XYZmm = voxel2mm(XYZ,V.mat);

	% for compatibility with SPM struct and cluster analysis
	% ------------------------------------------------------
	V.voxSize = [hdr.xsize hdr.ysize hdr.zsize]';
	V.u = crit_t;
	V.k = cl_size;
	V.title = V.fname;
	V.VOX = V.voxSize;		% compatible with VOL structure
    V.M = inV.mat;
    
	%eval(['save ' maskname '_struct V'])

return