function [chk, clusters] = check_spm_mat(mat1,mat2,clusters)
% mat1 is from clusters, mat2 is functional (imgs to extract)
%
% :Usage:
% ::
%
%     check_spm_mat(mat1,mat2,clusters)


chk = mat1 - mat2; chk = chk(:);
chk = chk(1:end-1);             % eliminate SPM scale factor
chk = any(chk);

if chk
    
    % we need to get the correct voxel coords from mat2 (funct) into
    % clusters, keeping the mm_coordinates the same!
    VOL.M = mat2;
    
    for i = 1:length(clusters)
        
        clusters(i).XYZ = mm2voxel(clusters(i).XYZmm,VOL)'; % functional img space, cluster mm coordinates
    
        clusters(i).Z = ones(1,size(clusters(i).XYZ,2));
        clusters(i).XYZmm = voxel2mm(clusters(i).XYZ,mat2);
        
        clusters(i).M = mat2;
        
        clusters(i).voxSize = diag(clusters(i).M(1:3,1:3)');
      
        clusters(i).numVox = size(clusters(i).XYZmm,2);
        
        % skip this, and clusters will have mm list of different length than
        % voxel list (XYZ).  data will be from voxel list.
        %SPM.XYZmm = voxel2mm(SPM.XYZ,VOL.M);    % slow, but gives unique voxels
    end
end

return
