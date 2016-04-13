function [indx, wh] = iimg_xyz2indx(xyz, xyztype, V)
% :Usage:
% ::
%
%     [indx, wh] = iimg_xyz2indx(xyz, input type:['mm' or 'vox'], [V: needed if 'mm'])

    if nargin > 1 && strcmp(xyztype, 'mm')
        xyz = mm2voxel(xyz, V);
    end

    indx = zeros(V.nvox, 1);
    wh = sub2ind(V.dim(1:3), xyz(:, 1), xyz(:, 2), xyz(:, 3));
    indx(wh) = 1;
end

