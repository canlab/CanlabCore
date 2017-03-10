function clusters = xyz2clusters(xyz,P)
% Converts a 3-column x, y, z list of mm coordinates to a clusters
% structure given P, the filename of an analyze .img file to provide
% dimensions and voxel sizes.
%
% :Usage:
% ::
%
%    function cl = xyz2clusters(xyz,P)
%
% Uses this info from the image:
%   - VOL.M     - spm-style mat matrix
%   - VOL.VOX   - voxel sizes
%
%   - SPM.Z     - now 1s; could stores values in the original image in clusters.Z
%
% The following is created internally:
%   - SPM.XYZmm - mm coords, you input these
%   - SPM.XYZ   - voxel coords
%

V = spm_vol(P);

%wh = strmatch('Thalamus',L3); whos wh


VOL.M = V.mat;
VOL.VOX = diag(V.mat(1:3,1:3));

fprintf(1,'Converting mm to voxels -> ');
XYZ = mm2voxel(xyz,VOL,2);  % unique reordered
fprintf(1,'%3.0f voxels -> ',size(XYZ,1));
fprintf(1,'Converting voxels to mm -> ');
%XYZ = unique(XYZ,'rows');
XYZ = XYZ';

SPM.XYZmm = voxel2mm(XYZ,VOL.M);
SPM.XYZ = XYZ;
SPM.Z = ones(1,size(XYZ,2));

fprintf(1,'Making clusters.\n');

%cl = tor_extract_rois([],SPM,VOL);
 
% stuff from tor_extract_rois -- omit where it gets each region separately
 
    % ----------------------------------------------------------------------------------
    % define each cluster as cell in array.
    % ----------------------------------------------------------------------------------
    clusters = [];
cl_index = spm_clusters(SPM.XYZ);
%cl_index = SPM.Z;    %:size(SPM.XYZ,2);
cl.title = 'xyz';

    for i = 1:max(cl_index)
      %  if verbose, fprintf(1,['\nExtracting cluster ' num2str(i) ' of ' num2str(max(cl_index))]),end
        a = find(cl_index == i);
        
        %if isfield(SPM,'u'), cl.title = SPM.title;,else, cl.title = 'Untitled';,end
        %if isfield(SPM,'u'), cl.threshold = SPM.u;,else, cl.threshold = NaN;,end
        cl(i).voxSize = VOL.VOX;
        cl(i).M = VOL.M;
        cl(i).name = [cl.title '_' num2str(i) '_' mat2str(size(a,2)) '_voxels'];
        cl(i).numVox = size(a,2);
        cl(i).Z = SPM.Z(a);
        cl(i).XYZmm = SPM.XYZmm(:,a);
        cl(i).XYZ = SPM.XYZ(:,a);
    
        cl(i).center = mean(cl(i).XYZ');
        if size(cl(i).XYZmm,2) > 1, cl(i).mm_center = center_of_mass(cl(i).XYZmm,cl(i).Z);  % mean(cl.XYZmm');
        else cl(i).mm_center = cl(i).XYZmm';
        end
        clusters = [clusters, cl];
        
    end
    
return
