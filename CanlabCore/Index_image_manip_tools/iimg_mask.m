function [dat, volInfo, masked_vol] = iimg_mask(maskindx,dat,varargin)
% Masks index image vector or filenames in dat
% With img vector or filenames in maskindx
% * Checks dimensions to make sure origins and vox. sizes are same
%
% :Usage:
% ::
%
%     [masked_indx, volInfo, masked_vol] = iimg_mask(maskindx,dat,[volInfo], [outname])
%
% Fastest way to mask an image:
% ::
%
%    [masked_indx_data, volInfo] = iimg_mask(maskindx,dat)
%
%    % or
%
%    [masked_indx, xyz, masked_vol] = iimg_mask(mask filename,image filenames)
%
% even slower: reconstruct a 3-D volume and return in masked_vol

 
[volInfo, maskindx] = iimg_read_img(maskindx,1);
[volInfod, dat] = iimg_read_img(dat);

% input volInfo overrides; in case maskindx is data vector instead of
% filename
if length(varargin) > 0 && ~isempty(varargin{1}), volInfo = varargin{1}; end


% deal with multiple images
% --------------------------------------------
[nvox,nimgs] = size(dat);


% reduce dat, if maskindx is reduced to in-mask only
% if both are image files, this accomplishes the masking
% --------------------------------------------
if nvox == volInfo.nvox && size(maskindx,1) ~= volInfo.nvox
    dat = dat(volInfo.wh_inmask,:);
end


% check sizes
% --------------------------------------------
if size(maskindx,1) - size(dat,1)
    error('Sizes of mask and image data do not match.');
end


if isstruct(volInfo) && isstruct(volInfod)
    % voxel sizes and origins
    m = [diag(volInfo.mat(1:3,1:3))' volInfo.mat(1:3,4)'];
    m2 = [diag(volInfod.mat(1:3,1:3))' volInfod.mat(1:3,4)'];
    
    if any(m-m2)
        error('Voxel sizes or origins for mask and image data do not match.');
    end
end
    
    
% masking
% --------------------------------------------
dat(~maskindx) = 0;


if nargout > 2
    % slower, reconstruct mask
    %if nargin < 3, error('You must enter volume structure V to reconstruct 3-D vol');, end
    
    for i = 1:nimgs
        if length(varargin) > 1
            outname = varargin{2};
            masked_vol(:,:,:,i) = iimg_reconstruct_3dvol(dat(:,i),volInfo, 'outname', outname);
            
        else
            masked_vol(:,:,:,i) = iimg_reconstruct_3dvol(dat(:,i),volInfo);
        end
        
    end
end

return
