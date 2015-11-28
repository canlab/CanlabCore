function [clusters,CLU,subclusters] = mask2clusters(P,varargin)
% Extracts clusters and con img data from mask
%
% Use with *mask_intersection.m*
%
% To get clusters but not extract data, enter only one argument.
%
% To get clusters and choose extraction imgs with the GUI, enter an empty [] 2nd argument.
%
% :Usage:
% ::
%
%    [clusters,CLU,subclusters] = mask2clusters(img mask file with voxels,[imgs to extract data from],[df])
%
%
% DOES *NOT* CONVERT BETWEEN DIFFERENT VOXEL SIZES AND POSITIONS BETWEEN IMNAMES AND SPM/VOL STRUCTS
%
% :See also: roi_probe
%
% If no imgs are entered, Z-scores are values from mask
%
% If df is entered, values in mask img are converted to Z-scores with spm_t2z.m
%
% If extract img names are empty and df is entered, assume we're using
% values from mask as t-values and convert to Z-scores
%
% WARNING: for spm2 compatibility, ABSOLUTE VALUES of voxel sizes are
% returned; e.g., ignores analyze flipping in SPM2.
%
% % Matlab 6.5/OSX bug gives seg fault or something if mask is too big.
%
% :Example:
% ::
%
%    cl = mask2clusters('myimage.img',[img string mtx],[]); % no z-score
%    conversion, extracts data from [img string mtx]
%
%    cl = mask2clusters('rob_tmap_0002_filt_t_3-05_k10_neg.img')
%
%    % This one works with already-loaded image data and a mat matrix:
%    V = spm_vol('rob_tmap_0002_filt_t_3-05_k10_neg.img'); dat = spm_read_vols(V);
%    cl = mask2clusters(dat,V.mat);
%
%
% ..
%    tor wager
%
%    modification 2/27/03
% ..

% ..
%    set up inputs
% ..
clusters = [];
df = []; imP = [];  resl = 0; vmat = [];


for i = 1:length(varargin)
    tmp = varargin{i};
    if isstr(tmp), imP = tmp;       % images to extract data from
    elseif isstruct(tmp), vmat = tmp.mat;  % it's an spm_vol V struct
    elseif isscalar(tmp), df = tmp; % df to convert from t to z-scores
    elseif all(size(tmp) > 1), vmat = tmp;      % it's a SPM .mat matrix to be used with raw data instead of string img name
    end
end

if isempty(P)
    P = spm_get(1,'*img','Select mask image with voxels to extract');
end

if ischar(P)
    P = deblank(P);
    
    if isempty(which(P)) && ~exist(P,'file')
        disp(['Cannot find file: ' P]);
        return
    end
    
    if ~isempty(which(P)), P = which(P); end
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% get coordinates and height values from the mask
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% img2voxel can take raw data, if vmat is entered; if string, vmat not
% used.
[CLU.XYZ,CLU.XYZmm,CLU.Z,CLU.V] = img2voxel(P,vmat);
CLU.XYZ = CLU.XYZ(1:3,:);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% change pseudo-T  or T values to Z scores
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if ~isempty(imP)
    if isempty(df) & length(varargin) > 1, df = size(imP,1) - 1;, end     % to get df
    if df == 0, df = [];, end
end

if ~isempty(df)
    disp(['Converting to Z scores based on ' num2str(df) ' df.'])
    [CLU.Z] = spm_t2z(CLU.Z,df);
else
    df = 0; 
    %disp('Saving values in mask file in clusters.Z (no z-score conversion)')
    %CLU.Z = ones(1,size(CLU.XYZ,2));
end

if size(CLU.Z,1) > size(CLU.Z,2), CLU.Z = CLU.Z';, end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% fill in other fields
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% adjust crit_t from 0 to 1 for compatibility with fixed_TSU (Talairach)
% needs to be a non-zero value to display on Talairach
% ------------------------------------------------------
CLU.crit_t = 1;

CLU.cl_size = 0;
CLU.df = df;

% for compatibility with SPM struct and cluster analysis
% ------------------------------------------------------
CLU.voxSize = diag(CLU.V.mat)'; 
CLU.voxSize = CLU.voxSize(1:3);
CLU.VOX = CLU.voxSize;		% compatible with VOL structure
CLU.M = CLU.V(1).mat;

CLU.u = CLU.crit_t;
CLU.k = CLU.cl_size;

if isfield(CLU.V,'fname')
    [a,b,c] = fileparts(CLU.V.fname);
    CLU.title = [b];
else
    CLU.title = 'Analyze image file.';
end

if resl, imP = spm_get(Inf,'*img','Select imgs for data extraction.');
else, % input imP              % to get P
end

%if ~exist(P(1,:))
%    warning('Image files cannot be found in ind. subject directories: Path has changed?')
    % tries to re-find P - only works for my specific directory 
%    for i = 1:size(P,1)
%        [d,f,e] = fileparts(P(1,:));
%        [dummy,d] = fileparts(d);
%        newP{i,1}=fullfile('..',d,[f e]);
%    end
%    P = cell2mat(newP);
%    disp(['Found: ' P(1,:) ' etc.'])
%end

if isempty(CLU),
    warning('EMPTY ... mask2clusters found no eligible voxels.')
    clusters = [];
    return
end

if isempty(CLU.XYZ) | isempty(CLU.Z)
    warning('EMPTY ... mask2clusters found no eligible voxels.')
    clusters = [];
    return
end

[clusters] = tor_extract_rois(imP,CLU,CLU,1);
for i = 1:length(clusters), 
    clusters(i).P = P;, 
    clusters(i).imP = imP;, 
    if size(imP,1) == 1 & df == 0,
        %disp('Saving values in mask file in clusters.Z')
        clusters(i).Z = clusters(i).all_data;
    end
    
    % for SPM2 compatibility
    clusters(i).voxSize = abs(clusters(i).voxSize);
end

if nargout > 2, [subclusters] = cluster_princomp(clusters); end

return
