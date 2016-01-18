function in = scn_mat_conform(in)
% :Usage:
% ::
%
%     function in = scn_mat_conform(in)
%
% sets flipping to 0 (no flip) in SPM2 and adjusts mat file accordingly
% input in spm-style mat file or struct with .mat or .M fields
%
% ..
%    tor wager, dec 06
% ..

    global defaults

    if isstruct(in)
        if isfield(in,'mat')
            in.mat = scn_mat_conform(in.mat);
        end

        if isfield(in,'M')
            in.M = scn_mat_conform(in.M);
        end

        return
    end



    if in(1) < 0
        disp('Warning: Image has negative x voxel size, indicating ''flipped'' in SPM.')
        disp('This will be changed to positive. This program does not do any image flipping.')

        in(1) = abs(in(1));     % voxel size
        in(1,4) = -(in(1,4));   % origin offset
    end

    if isempty(defaults)
        spm_defaults
    end

    if defaults.analyze.flip
        disp('Warning: Setting defaults.analyze.flip to 0.  No flipping.')
        defaults.analyze.flip = 0;
    end


    return
