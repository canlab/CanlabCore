function [mask, mask2] = clusters2mask2011(cl, varargin)
% Returns 3-D mask of voxels included in a cl structure.
%
% :Usage:
% ::
%
%     [mask, mask2] = clusters2mask2011(cl, [dim])
%
% Mask values are coded with the cluster index. That is, the non-zero entries
% in mask are integers that reflect the index of the unique contiguous cluster to
% which each voxel belongs.
%
% Any non-zero value indicates membership in a cluster
% Uses VOXEL values from cl, so define cl in the space you wish to have for
% mask first!
%
% :Optional:
%
%   A 2nd argument will be treated as 'dim', dimensions in voxels of the
%   new mask image.  If empty, uses max value in cluster to determine
%   automatically, but then the mask may not match image dimensions desired
%   for .img/.nii reading/writing purposes.
%
%   If a second output is requested, the second output (maskz) is a mask
%   like the first, but the values in the mask reflect the numeric value
%   stored in the cl.Z field (whether Z-scores or other values, depending on
%   how cl is constructed.)
%
% ..
%    Copyright 2011 Tor Wager
% ..


mask = [];
mask2 = [];

if isempty(cl), return, end

if isa(cl, 'region') && ~any(cat(1, cl.numVox))
    [mask, mask2] = deal(zeros(cl(1).dim));
    return
end

if ~isfield(cl(1), 'M') || isempty(cl(1).M)
    %error('cl(1).M must have spm-style mat info, i.e., as returned in spm_vol.m V.mat');
    % do nothing; use voxel values alone
    % if we have .mat file entered from some other source, could use that
    % to reconstruct.  but probably better to recon in voxel dims anyway,
    % and then resample to the grid defined by the new .mat information.
    % Then we will have voxels filled in appropriately when upsampling.
    % So the .M field is really not necessary at all.
end

XYZ = cat(2, cl(:).XYZ);

if length(varargin) == 0
    % get dim automatically
    dim = max(XYZ, [], 2);
    
else
    % we have dim entered
    dim = varargin{1};
    
    if any(size(dim) - [3 1])
        dim = dim';
    end
    
    if any(size(dim) - [3 1])
        error('dim argument must be 3 x 1 vector of image dimensions, in voxels');
    end
      
    
end


for i = 1:length(cl)
    wh_cl{i} = repmat(i, size(cl(i).XYZ, 2), 1);
end
wh_cl = cat(1, wh_cl{:});


ind = sub2ind(dim', XYZ(1, :), XYZ(2, :), XYZ(3, :));

[mask, mask2] = deal(zeros(dim'));
mask(ind) = wh_cl;

if nargout > 1
    mask2(ind) = cat(2, cl(:).Z);
end

end

