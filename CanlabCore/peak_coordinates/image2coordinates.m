function XYZmm = image2coordinates(img)
% Threshold a statistic image and turn it into a list of x, y, z
% coordinates (in mm "world space") with SPM
%
% :Usage:
% ::
%
%     XYZmm = image2coordinates(img)
%
% :Examples:
% ::
%
%    img = 'h25_aerger.img';
%    XYZmm = image2coordinates(img)
%    figure;
%    plot3(XYZmm{1}(1, :)', XYZmm{1}(2, :)', XYZmm{1}(3,:)', 'ko');
%    addbrain
%    axis image
%
% Example with multiple images:
% ::
%
%    % list all images in dir
%    imgs = filenames(fullfile(pwd, '*img'), 'char', 'absolute')
%    XYZmm = image2coordinates(imgs); 
%    % save images names and coordinates
%    % imgs2 has cell array of image names, without full paths
%    imgs2 = filenames(fullfile(pwd, '*img')); 
%    save silke_senders_xyz_coordinates imgs imgs2 XYZmm

dat = fmri_data(img);

% arbitrary threshold of 1 or greater, based on histogram
%dat = threshold(dat, [1 Inf], 'raw-between');

% threshold is xth (e.g., 99th) percentile across all the images
thr = prctile(double(dat.dat(:)), 99);
dat = threshold(dat, [thr Inf], 'raw-between');

dat = replace_empty(dat);

n = size(dat.dat, 2);  % number of images

XYZmm = cell(1, n);

for i = 1:n
    
    % values in image(s) - need for spm_max
    X = double(dat.dat(:, i));
    wh = X ~= 0;
    
    % remove empty voxels
    X(X == 0) = [];
    
    % skip if empty (no voxels)
    if isempty(X), continue, end
    
    % otherwise
    xyz = dat.volInfo.xyzlist(wh, :)';
    
    % core SPM function.  M is maxima in voxel coords
    [N Z M] = spm_max(X, xyz);
    
    % put into mm:
    XYZmm{i} = voxel2mm(M, dat.volInfo.mat);
    
end

end % function
