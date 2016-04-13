function slice_data = get_ax_slice(imgs, slice_num)
% Get an axial slice
%
% :Usage:
% ::
%
%     function slice_data = get_ax_slice(imgs, slice_num)
%
% :Inputs:
%
%   **imgs:**
%        img filenames or spm_vols eg. 'spmT_0004.img';
%
%   **slice_num:**
%        slice number eg. 31
%
% :Output:
%
%   **slice_data:**
%        unprocessed slice data
%

    global defaults
    
    
    V = spm_vol(imgs);

    % To figure out the appropriate C matrix to use in spm_slice_vol
    % use the fact that
    %
    % [x,y,z,1]' = C [x',y',0,1]
    %
    % where [x',y'] define coordinates in your plane of interest.
    % A voxel in the volume z=[x,y,z,1]' which corresponds to
    % a pixel in this slice is then
    % a linear combination of [x',y']. For example, to
    % get a transverse slice the z-value is always fixed at p.
    % So the 3rd row of C is [0 0 1 p].

    % Transverse slice
    %==========================================================
    C = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
    DIM = V(1).dim(1:2);

    C(3,4) = slice_num;
    %C(3,4)=-p;

    % img = rot90(spm_slice_vol(V,C,DIM,0));
    % img = spm_slice_vol(V,inv(C),DIM,0);
    for i=1:length(V)
        slice_data(:,:,i) = spm_slice_vol(V(i), C, DIM, 0);
    end

    slice_data = squeeze(slice_data);
end
