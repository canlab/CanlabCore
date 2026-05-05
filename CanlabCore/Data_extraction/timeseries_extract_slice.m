function sl = timeseries_extract_slice(V, sliceno, orientation)
% For a given set of image names or memory mapped volumes (V)
% extracts data from slice # sliceno and returns an X x Y x time
% matrix of data.
%
% :Usage:
% ::
%
%     function sl = timeseries_extract_slice(V, sliceno, orientation)
%
% :Inputs:
%
%   **V:**
%        image filenames (char) or spm_vol memory-mapped volumes
%
%   **sliceno:**
%        slice number (voxel index) to extract
%
%   **orientation:**
%        'axial' (default), 'sagittal', or 'coronal'
%
% :Output:
%
%   **sl:**
%        X x Y x time matrix of slice data
%

    if ischar(V)
        V = spm_vol(V);
    end

    if(~exist('orientation', 'var') || isempty(orientation))
        orientation = 'axial';
    end

    switch(orientation)
        case 'axial'
            % Map output pixel (x,y) -> voxel (x, y, sliceno)
            mat = spm_matrix([0 0 sliceno]);
            for i = 1:length(V)
                sl(:,:,i) = spm_slice_vol(V(i), mat, V(i).dim(1:2), 0);
            end

        case 'sagittal'
            % Map output pixel (x,y) -> voxel (sliceno, x, y)
            % Matrix maps: x_vox=sliceno (const), y_vox=x_pix, z_vox=y_pix
            mat = [0 0 0 sliceno; 1 0 0 0; 0 1 0 0; 0 0 0 1];
            for i = 1:length(V)
                sl(:,:,i) = spm_slice_vol(V(i), mat, V(i).dim(2:3), 0);
            end

        case 'coronal'
            % Map output pixel (x,y) -> voxel (x, sliceno, y)
            % Matrix maps: x_vox=x_pix, y_vox=sliceno (const), z_vox=y_pix
            mat = [1 0 0 0; 0 0 0 sliceno; 0 1 0 0; 0 0 0 1];
            for i = 1:length(V)
                sl(:,:,i) = spm_slice_vol(V(i), mat, [V(i).dim(1) V(i).dim(3)], 0);
            end

        otherwise
            error('Unknown orientation: %s\n', orientation);
    end
end
