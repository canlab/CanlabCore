% :Usage:
% ::
%
%     [g, gslice, stdslice, s]] = tor_global(P, [maskname])
%
% Compute global image values within an optional mask
% alternative to spm_global
%
% Pass in a file name for mask:
% ::
%
%    g = tor_global(P, 'my_gray_mask.img')
%
% Use Worsley's fmri_mask_thresh to implicitly mask based on first image:
% ::
%
%    g = tor_global(P, 1);
%
% No masking, use whole image
% ::
%
%    g = tor_global(P)
%
% :Outputs:
%
%   **g:**
%         array mean intensity for each image
%
%   **gslice:**
%         2D array with mean intensity for each slice
%
%   **stdslice:**
%         2D array with std for each slice
%
%   **s:**
%         array std for each image

function [g, gslice, stdslice, s] = tor_global(P, varargin)
    
    switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

        case {'SPM5', 'SPM8', 'SPM12'}
            % spm_defaults is a function
            spm_defaults()

        otherwise
            % unknown SPM
            disp('Unknown version of SPM!');
            spm_defaults()
    end

    fprintf('getting globals.\t');
    
    V = spm_vol(P);
    num_vols = length(V);
    mask = [];

    if ~isempty(varargin)
        % we have some mask info
        % If mask is left empty, then we will just not mask at all.
        maskinfo = varargin{1};

        if isscalar(maskinfo) && maskinfo  % It's just a 1 or something, implicit mask
            % ----------------------------------------------------------------------------------
            % * create mask from FIRST IMAGE and mask all subjects'  values
            % so that we show only those voxels existing for all subjects.
            % ----------------------------------------------------------------------------------
            % it's a threshold
            if ~exist('fmri_mask_thresh.m', 'file')
                error('You must have fmri_mask_thresh.m from Keith Worsley''s FMRISTAT package.');
            end
            mask_thresh = fmri_mask_thresh(V(1).fname);
            mask_dat = spm_read_vols(V(1));
            mask = mask_dat > mask_thresh;  % create mask data
            clear mask_dat

        elseif ischar(maskinfo)
            % ----------------------------------------------------------------------------------
            % * load mask, from FILE and mask all subjects'  values
            % so that we show only those voxels existing for all subjects.
            % ----------------------------------------------------------------------------------          %

            fprintf('masking volumes.\t')
            % maskinfo is a filename
            mask = scn_map_image(maskinfo, V(1).fname);
            mask = double(mask ~= 0 & ~isnan(mask));
        end

    end

    g = zeros(num_vols, 1);
    gslice = [];


    new_status = sprintf('Image %4d', 0);
    fprintf(new_status);
    for i = 1:num_vols
        old_status = new_status;
        new_status = sprintf('Image %d/%d', i, num_vols);
        erase_and_display(old_status, new_status);

        [g(i), gslice(:,i), stdslice(:, i), s(i)] = get_global(V(i), mask);
    end

    fprintf('\nFinished computing global on %d images.\n', num_vols);
end



% V is a 3D image
% returns:
%   g: mean intensity for the image
%   gslice: array with mean intensity for each slice,
%   stdslice: array with std for each slice
%   s: std for image
function [g, gslice, stdslice,s] = get_global(V, mask)
    dat = spm_read_vols(V);

    if ~isempty(mask), dat = dat .* mask; end


    % slice-by-slice globals
    z = size(dat, 3); %num of slices
    gslice = zeros(z, 1);

    for i = 1:z
        slice = dat(:,:,i);
        slice = slice(:);
        slice(isnan(slice) | slice == 0) = [];

        if ~isempty(slice)
            gslice(i) = mean(slice); %mean value for slice
            stdslice(i) = std(slice); % STD for slice
        else
            gslice(i) = 0;
            stdslice(i) = 0;
        end
    end

    dat = dat(:);
    dat(isnan(dat) | dat == 0) = [];

    g = mean(dat); %mean value for image
    
    s = std(dat); %STD for image
    
end



