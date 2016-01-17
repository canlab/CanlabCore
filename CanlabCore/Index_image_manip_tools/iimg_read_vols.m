function [Y, XYZ] = iimg_read_vols(V, mask)
% Drop-in replacement for spm_read_vols
% Primary difference is that if the images have the same voxel size but are NOT resliced
% (e.g., have differing affine matrices), this will handle reading each individually and
% putting them together.


    if nargin<2, mask = 0; end
    if nargin<1, error('insufficient arguments'); end

    if length(V)>1 && any(any(diff(cat(1, V.dim), 1, 1), 1))
        error('images don''t all have the same dimensions');
    end
    if any(any(any(diff(cat(3, V.mat), 1, 3), 3)))
        %error('images don''t all have same orientation & voxel size');
        for i = 1:length(V)
            Y(:,:,:,i) = spm_read_vols(V(i), mask);
        end
        
        if nargout > 1
            [R, C, P] = ndgrid(1:V(1).dim(1), 1:V(1).dim(2), 1:V(1).dim(3));
            RCP = [R(:)';C(:)';P(:)'];
            clear R C P
            RCP(4,:) = 1;
            XYZ = V(1).mat(1:3,:)*RCP;
        end
    else
        [Y, XYZ] = spm_read_vols(V, mask);
    end
end
