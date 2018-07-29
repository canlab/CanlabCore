function obj = reparse_contiguous(obj, varargin)
% Re-construct list of contiguous voxels in an image based on in-image
% voxel coordinates.  Coordinates are taken from obj.volInfo.xyzlist.
% Results are saved in obj.volInfo.cluster.
% xyzlist can be generated from iimg_read_img, and is done automatically by
% object-oriented fMRI image classes (fmri_image, image_vector,
% statistic_image)
%
% :Usage:
% ::
%
%    obj = reparse_contiguous(obj, ['nonempty'])
%
% If 'nonempty' is entered as an optional argument, will use only voxels
% that are non-zero, non-nan in the first column of obj.dat.
%
% The statistic_image object version of reparse_contiguous uses 
% the significance of the first image in the object (obj.sig(:, 1)) as a
% filter as well, so clustering will be based on the latest threshold applied.
% it is not usually necessary to enter 'nonempty'.
%
% :Example:
% ::
%
%    % Given timg, a statistic_image object:
%    test = reparse_contiguous(timg, 'nonempty');
%    cl = region(test, 'contiguous_regions');
%    cluster_orthviews(cl, 'unique')
%
% ..
%    Copyright tor wager, 2011

%    Programmers' notes:
%    6/22/14: Tor changed behavior to use .sig field. Returns clusters of
%    in-mask, significant (.sig) voxels only.
% %  7/2018 - tor - fixed bug with statistic_image handling in some cases
%             tor - clusters not in .sig were not zeroed out - now they are
% ..

wh = true(size(obj.volInfo.cluster));   %obj.volInfo.wh_inmask;
obj.volInfo(1).cluster = zeros(size(wh));

obj = replace_empty(obj);

% restrict to voxels with actual data if desired
if any(strcmp(varargin, 'nonempty'))
        
    wh = all(obj.dat ~= 0, 2) & all(~isnan(obj.dat), 2);
    
end

% 6/22/13 Tor Added to enforce consistency in objects across usage cases
if isempty(obj.sig) || (numel(obj.sig) == 1 && ~obj.sig)
    obj.sig = true(size(obj.dat));
end

% .cluster can be either the size of a reduced, in-mask dataset after removing empties
% or the size of the full in-mask dataset that defined the image
% (volInfo.nvox).  We have to switch behavior according to which it is.
if size(obj.volInfo(1).cluster, 1) == size(obj.volInfo(1).xyzlist, 1)

    wh = logical(obj.sig(:, 1));
    
    newcl = spm_clusters(obj.volInfo(1).xyzlist(wh, :)')'; % 6/22/14 tor changed to use .sig
    
    obj.volInfo(1).cluster(wh) = newcl;                     % tor changed to use .sig
    
    obj.volInfo(1).cluster(~wh) = 0;                        % zero to clusters outside .sig mask
    
else % full in-mask -- assign to correct voxels
    
    wh = wh & logical(obj.sig(:, 1));

    obj.volInfo(1).cluster(wh) = spm_clusters(obj.volInfo(1).xyzlist(wh, :)')';
    
     obj.volInfo(1).cluster(~wh) = 0;                       % zero to clusters with no data or outside .sig mask
end

obj = remove_empty(obj);

end



