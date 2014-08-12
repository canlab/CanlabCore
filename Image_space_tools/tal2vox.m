function vox=tal2vox(tal,VOL)
% converts from talairach coordinate to voxel coordinate
% based on variables from SPM.M (passed here for 
% faster operation)
% e.g., foo=tal2vox([-30 28 -30], VOL)

if(isfield(VOL, 'M'))
    M = VOL.M;
elseif(isfield(VOL, 'mat'))
    M = VOL.mat;
else
    error('Error in %s: VOL does not have an "M" or a "mat" field.', mfilename);
end

vox=[0 0 0];
vox(1)=(tal(1)-M(1,4))/M(1,1);
vox(2)=(tal(2)-M(2,4))/M(2,2);
vox(3)=(tal(3)-M(3,4))/M(3,3);

return

