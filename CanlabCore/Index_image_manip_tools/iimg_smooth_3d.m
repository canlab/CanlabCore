function w = iimg_smooth_3d(w, volInfo, sfwhm, varargin)
% Smooth 3-D images stored in columns with FWHM in mm of sfwhm
%
% :Usage:
% ::
%
%     function w = smooth_3d(w, volInfo, sfwhm, [badvox])
%
% NOTE: SMOOTHING KERNEL MAY BE IN VOX, AS VOL INFO IS NOT PASSED IN
%
% :Inputs:
%
%   **w:**
%        Take v x n matrix w, and smooth columns (images), returning
%        v x n matrix again
%
%   **volInfo:**
%        is volume info structure from iimg_read_img.m
%
% :Optional: If v is smaller than the original image because some voxels
% were removed, enter logical vector badvox, and the missing voxels
% will be filled with zeros.
%
% NOTE: 4-D version is horribly slow and memory intensive:
%
% :Example:
% ::
%
%    wvol = iimg_reconstruct_vols(w', fmri_data_obj.mask.volInfo);

badvox = 0;
if ~isempty(varargin), badvox = varargin{1}; end

% transpose, and ...
% insert zeros back into bad vox
if any(badvox)
    w = zeroinsert(badvox, w);
end

n = size(w, 2);

try
    % reconstruct each image separately and smooth it,
    % then convert back.
    
    for i = 1:n
        fprintf('\b\b\b\b%04d', i);
        
        wvol = iimg_reconstruct_vols(w(:, i), volInfo);
        
        spm_smooth(wvol, wvol, [sfwhm sfwhm sfwhm]);
        
        % go back to vector
        wvec = wvol(volInfo.wh_inmask);
        
        w(:, i) = wvec;
        
    end
    
catch
    disp('Error with volume reconstruction. Mask in fmri_data obj not properly defined?');
    disp('Stopped in debugger so you can check...');
    keyboard
end

if any(badvox)
    
    w(badvox, :) = [];
    
end

end % function

