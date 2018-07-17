function obj = reparse_contiguous(obj, varargin)
% Re-construct list of contiguous voxels in an image based on in-image
% voxel coordinates.  Coordinates are taken from obj.volInfo.xyzlist.
%
% Results are saved in obj.volInfo.cluster.
%
% xyzlist can be generated from iimg_read_img, and is done automatically by
% object-oriented fMRI image classes (fmri_image, image_vector,
% statistic_image)
%
% If 'nonempty' is entered as an optional argument, will use only voxels
% that are non-zero, non-nan in all columns of obj.dat.
%
% :Usage:
% ::
%
%    obj = reparse_contiguous(obj, ['nonempty'])
%
% ..
%    Copyright tor wager, 2011
%
%    Programmers' notes:
%
%    Edited 1/27/13 by tor to use all columns when calculating 'nonempty'
%                also fixed bug - was not using 'nonempty' input in some
%                cases
% ..
%
% .cluster and .xyzlist should both always be length v in-mask voxels
% if 'nonempty' is entered, then .dat should be length v in-mask voxels too
%
% NOTES: this will only work if xyzlist in volInfo has only voxels for
% contiguous clusters. If the mask in volInfo contains the whole brain,
% use 'nonempty' to exclude empty data values when
% redefining clusters.

wh = true(size(obj.volInfo.cluster));   %obj.volInfo.wh_inmask;

if isempty(wh)
    
    disp('Warning: image_vector.reparse_contiguous: obj.volInfo.cluster is empty. Bad object...rebuilding.');

    % No clusters - rebuild cluster index
    obj.volInfo.cluster = spm_clusters(obj.volInfo.xyzlist')';
    
    
end

if isempty(obj.dat)
    
    disp('Warning: image_vector.reparse_contiguous: obj.dat is empty. Bad object...');
    disp('Stopping in debugger so you can check object.');
    keyboard
    
end



% restrict to voxels with actual data if desired
if any(strcmp(varargin, 'nonempty'))
    
    obj = replace_empty(obj,'voxels');  % changed to 'voxels' only. SG 2/23/18
    wh = all(obj.dat ~= 0, 2) & all(~isnan(obj.dat), 2); % this will expand wh to full in-mask dataset, .nvox
    
end

obj.volInfo(1).cluster = zeros(size(wh));


% .cluster can be either the size of a reduced, in-mask dataset after removing empties
% or the size of the full in-mask dataset that defined the image
% (volInfo.nvox).  We have to switch behavior according to which it is.

if size(obj.volInfo(1).cluster, 1) == size(obj.volInfo(1).xyzlist, 1)
    % if we have in-mask voxel list only in .cluster, .xyzlist
    
    if length(wh) ~= length(obj.volInfo(1).cluster)
        error('Bug: length(wh) does not match number of elements in volInfo.cluster.');
    end
    
    newcl = spm_clusters(obj.volInfo(1).xyzlist(wh, :)')';
    obj.volInfo(1).cluster = zeros(size(wh));
    obj.volInfo(1).cluster(wh) = newcl;
    
elseif size(obj.volInfo(1).xyzlist, 1) ~= length(wh)
    
    disp('obj.volInfo.xyzlist voxels:')
    size(obj.volInfo(1).xyzlist, 1)
    
    disp('Voxels in image list to re-parse:')
    length(wh)
    
    error('You have created an invalid image object somewhere along the line.');
    
else   % full in-mask -- assign to correct voxels
    
    obj.volInfo(1).cluster(wh) = spm_clusters(obj.volInfo(1).xyzlist(wh, :)')';
    
end

end % function




