function mask = voxel2mask(voxels,maskdims)
% function mask = voxel2mask(voxels, x y z mask dimensions)
% voxels:
% 3 column vectors
% [i j k] = row, column, slice
% [x y z] in brain if brain is in analyze format
% (x is rows, y is columns, z is slices)
% Tor Wager, 10/17/01

[n, m] = size(voxels);
wh_bad = false(n, 1);

% Check for illegal voxels
if m ~= 3, warning('voxel2mask: Illegal voxel input list', 'Voxels must be k x 3 matrix'); end

mv = max(voxels);
if any(mv > maskdims)
    warning('voxel2mask: Illegal voxel input list', 'Voxels outside mask, Voxel indices > mask dims');
    wh_bad = any(voxels > mv(ones(n, 1), :), 2);
end

mv = min(voxels);
if any(mv < 1)
    warning('voxel2mask: Illegal voxel input list', 'Voxels outside mask, Voxel indices < 1');
    wh_bad = [wh_bad | any(voxels < ones(n, 3), 2)];
end

if ~isempty(wh_bad)
    wh_bad = find(wh_bad);
    disp('You need to check your images!');
    disp('Offending voxel numbers: '); disp(wh_bad)
    disp('Offending voxel coordinates: '); disp(voxels(wh_bad, :))
    voxels(wh_bad, :) = [];
end

mask = zeros(maskdims);

for i = 1:size(voxels, 1)
    mask(voxels(i,1),voxels(i,2),voxels(i,3)) = 1;
end

mask = double(mask);

return