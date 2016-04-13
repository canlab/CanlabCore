function dat = mask_image(img, mask, outname, varargin)
% :Usage:
% ::
%
%     masked_dat = mask_image(img, mask, outname, ['reverse'])
%
% Mask an image file (img) with a mask (mask), and save in outname.
% zero and NaN values are considered invalid values, and so voxels with
% these values are considered "excluded"
%
% :Optional Inputs:
%
%   **'minmask':**
%        mask values less than next argument are excluded
%
%   **'maxmask':**
%        mask values greater than next argument are excluded
%
%   **'minimg':**
%        img values less than next argument are excluded
%
%   **'maximg':**
%        img values greater than next argument are excluded
%
%   **'abs':**
%        impose min/max thresholds based on absolute values
%
%   **'reverse':**
%        make "reverse mask," including only previously excluded
%        areas (values of zero or NaN)
%
%        Note: applies to mask, not img values
%        So values with 0 in img will always be 0, whether
%       standard or "reverse" mask is used.
%
%   **'binary':**
%        make the image values binary (i.e., create a new mask)
%
% This function can handle images of different dimensions.  The output
% image will use the dimensions of img.
%
% :Examples:
% ::
%
%    % Create an image with non-zero numbers only where p-values in an image are greater than .05
%    img = 'X-M_pvals.img'
%    mask = 'X-M_pvals.img';
%    maxmask = .05
%    outname = 'notX_p05.img';
%    mask_image(img, mask, outname, 'reverse', 'maxmask', maxmask);
%    spm_image('init', outname);
%
%    mask_image(my_mean_image, 'functional_mask.img', ...
%                  'functional_mask.img', 'minimg', cutoff, 'abs');
%
%
%    mask_image('n15_avgpet.img',EXPT.mask,'n15_avgpet_brain.img');
%

    doreverse = 0;
    minmask = 0;
    maxmask = NaN;
    minimg = 0;
    maximg = NaN;
    doabs = 0;
    dobinary = 0;
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'reverse', doreverse = 1;

                case {'abs', 'absolute'}, doabs = 1;
                    
                case {'minmask'}, minmask = varargin{i + 1};

                case {'maxmask'}, maxmask = varargin{i + 1};

                case {'minimg'}, minimg = varargin{i + 1};

                case {'maximg'}, maximg = varargin{i + 1};

                case {'binary'}, dobinary = 1;
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end


    V = spm_vol(img); v = spm_read_vols(V);

    %M = spm_vol(mask); %m = spm_read_vols(M);
    m = scn_map_image(mask, img);

    % apply threshold(s), if any
    % -------------------------------------------
    if doabs, absstr = ' (absolute value of) '; else absstr = [' ']; end
    
    if minmask
        fprintf('Excluding voxels with%smask values less than %3.4f from mask area.\n', absstr, minmask);

        if doabs, to_use = abs(m); else to_use = m; end 
        m(to_use < minmask) = 0;
    end

    if ~isnan(maxmask)
        fprintf('Excluding voxels with%smask values greater than %3.4f from mask area.\n', absstr, maxmask);

        if doabs, to_use = abs(m); else to_use = m; end 
        m(to_use > maxmask) = 0;
    end

    if minimg
        fprintf('Excluding voxels with%simg values less than %3.4f from mask area.\n', absstr, minimg);

        if doabs, to_use = abs(v); else to_use = v; end 
        v(to_use < minimg) = 0;
    end

    if ~isnan(maximg)
        fprintf('Excluding voxels with%simg values greater than %3.4f from mask area.\n', absstr, maximg);

        if doabs, to_use = abs(v); else to_use = v; end         
        v(to_use > maximg) = 0;
    end



    % Mask
    % -------------------------------------------
    if doreverse
        % Reverse mask
        fprintf('Creating reverse mask.\n');
        dat = double(v .* (abs(m) == 0 | isnan(m)));
    else
        % Standard mask
        dat = double(v .* (abs(m) > 0));

    end
    V.fname = outname;

    if dobinary
        disp('Making output image binary')
        dat = double(abs(dat) > 0) & ~isnan(dat);
    end
    
    fprintf('Writing %s\n', V.fname);

    spm_write_vol(V,dat);

    % end

    return

