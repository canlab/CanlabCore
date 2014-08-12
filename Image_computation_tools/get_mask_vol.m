function mask = get_mask_vol(Pf,Pm)
% function mask = get_mask_vol(Pf,Pm)
% 
% Tor Wager
% Returns 3-D volume of binary mask values,
% in the space of the functional image Pf.
%
% Pf    filename of functional image
% Pm    filename of mask image (should be 1's and 0's)

    Vm = spm_vol(Pm);
    Vf = spm_vol(Pf);
    if diag(Vm.mat(1:3,1:3)) ~= diag(Vf.mat(1:3,1:3))   % unequal voxel sizes
        disp('Reslicing mask image to dims of functional img')
        P = reslice_imgs(Pf,Pm,0);
        [d f e] = fileparts(Pm);
        Pm = fullfile(d,['r' f e]);
        Vm = spm_vol(Pm);
    end
    mask = spm_read_vols(Vm);
    mask = double(mask > .5);

return

