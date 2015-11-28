function [m,V,cl] = clusters2mask(cl,V,varargin)
% This function has 2 modes!  If V is a structure:
%
% :Usage:
% ::
%
%    [m,V,cl] = clusters2mask(cl,V,[opt: write Z-scores])
%
% Converts clusters structure to a mask image, given V structure with V.mat
% field.  V.mat is an SPM mat file. V.dim is dims of image uses cl.XYZmm
% m is mask img data, V is mask vol info.
%
% Also replaces cl.XYZ (voxels)
%
% If V is a vector of mask dimensions:
% converts clusters to mask image using existing XYZ and dims of mask
%
% :See also: voxels2mask, for a faster function that uses XYZ voxel coords
%    
% :Example:
% ::
%
%    % Save an image file with just one cluster from a set (#7 in this ex.)
%    cl = mask2clusters('roi_group1.img');
%    V = spm_vol('roi_group1.img');  % we need .mat and .dim from this, or
%
%    % just dim
%    [m,V,cl] = clusters2mask(cl(7),struct('mat',cl(1).M,'dim',V.dim));
%    %or
%    [m,V,cl] = clusters2mask(cl(7),V);
%
%    clusters2mask(cl,struct('mat',V.mat,'dim',V.dim),0,'spm2_hy.img');
%
%    %for SPM5:
%    clusters2mask(cl,struct('mat',V.mat,'dim',V.dim, 'dt', V.dt),0,'spm2_hy.img');
%    clusters2mask(cl,
%    struct('mat',MC_Setup.volInfo.mat,'dim',MC_Setup.volInfo.dim, 'dt', MC_Setup.volInfo.dt),0,'acc_roi_mask.img');
%
% ..
%    tor wager
% ..

    global defaults
    
    fname = 'clustermask.img';
    if length(varargin) > 1, fname = varargin{2}; end

    doZ = 0;
    if length(varargin) > 0, doZ = varargin{1}; end

    if isstruct(V)
        %%% MM way

        m = zeros(V.dim(1:3));

        if isfield(V,'mat')
            %mymat = V.mat; do nothing
        elseif isfield(V,'M');
            V.mat = V.M;
        else
            error('Structure must have a mat or M field.');
        end

        for i = 1:length(cl)

            cl(i).XYZ = mm2voxel(cl(i).XYZmm,struct('M',V.mat),1)';

            % eliminate out of range


            ind = sub2ind(V.dim(1:3),cl(i).XYZ(1,:)',cl(i).XYZ(2,:)',cl(i).XYZ(3,:)');

            if doZ, m(ind) = cl(i).Z;,else, m(ind) = 1; end

        end

        switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

        case 'SPM5'
            % spm_defaults is a function
            
            if(isempty(defaults))
                spm_defaults();
            end
            
            if ~isfield(V, 'dt')
                warning('Using default datatype; Enter dt field in input structure to use yours.');
                V.dt(1) = spm_type('float32');
                V.dt(2) = 1;
            end
        end
    
        V.fname = fname;
        spm_write_vol(V,double(m));
        disp(['Written ' fname])
        cl = mask2clusters(fname);

    else
        %%% VOXEL way
        for i = 1:length(cl)
            mask(:,:,:,i) = voxel2mask(cl(i).XYZ',V);
        end
        m = sum(mask,4);
        m = double(m > 0);
    end

    return
