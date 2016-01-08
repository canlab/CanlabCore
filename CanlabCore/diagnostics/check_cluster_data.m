function check_cluster_data(cl)
% loads the first 5 images from the first voxel
%
% :Usage:
% ::
%
%     check_cluster_data(cl)
%
% ..
%    tor wager
%    needs Tom's orthviews
% ..


%cl(1).imnames(1:5,:)
cl(1).XYZ(:,1:5), cl(1).XYZmm(:,1:5),
spm_orthviews('Reposition',cl(1).XYZmm(:,1))
%cl(1).all_data(1:5,1)
cl(1).raw_data(1:5,1,2)

return
