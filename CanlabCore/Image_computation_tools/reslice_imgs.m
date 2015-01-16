% function [P, reslicedImgs] = reslice_imgs(sampleTo, resliceThis, [domask], [overwrite])
%
% arguments are file names of .img files
% if empty, select from GUI
%
% if domask, recalculates a 1 or 0 mask for each image file in resliceThis
% if overwrite, overwrite the original instead of prepending an 'r'
%
% Example:
% Reslice a mask image into the space of some functional images, and move
% to the current directory
% -----------------------------------------------------------------------
% [tmp, maskname] = reslice_imgs(image_names(1, :), maskname, 1);
% eval(['!mv ' maskname ' ./'])
% eval(['!mv ' maskname(1:end-4) '.hdr ./'])

function [P, reslicedImgs] = reslice_imgs(sampleTo, resliceThis, varargin)
    domask = 0;
    overwrite = 0;

    % using b-spline or fourier really screwed up some of my masks, so trilinear may
    % be better!
    flags = struct('interp', 1, ... % b-spline
        'mask', 0, ...              % do not mask
        'mean', 0, ...              % do not write mean image
        'hold', -1, ...             % i don't think this is used anymore
        'which', 1, ...             % reslice 2nd-nth only
        'wrap', [0 0 0]' ...        % the default; don't know what this is
        );

    if length(varargin) > 0
        domask = varargin{1};
    end
    if length(varargin) > 1
        overwrite = varargin{2};
    end

    if isempty(sampleTo)
        sampleTo = spm_get(1, '*.img', 'Select image with desired dimensions', pwd, 0);
    end

    if isempty(resliceThis)
        resliceThis = spm_get(Inf, '*.img', 'Select image(s) to resample', pwd, 0);
    end

    P = str2mat(sampleTo, resliceThis);
    spm_reslice(P, flags)


    if domask
        disp('Re-thresholding masks...')

        for i = 1:size(resliceThis, 1)
            [d, f, e] = fileparts(deblank(resliceThis(i, :)));
            resliced_img = new_resliced_name(overwrite, d, f, e);
            spm_imcalc_ui(resliced_img, resliced_img, 'i1>0');
        end
    end

    for i = 1:size(P, 1)-1
        [d, f, e] = fileparts(deblank(P(i+1, :)));
        resliced_img = new_resliced_name(overwrite, d, f, e);
        
        if(overwrite)
            movefile(fullfile(d, ['r' f e]), resliced_img);
            movefile(fullfile(d, ['r' f '.hdr']), [resliced_img(1:end-4) '.hdr']);
            if(exist(fullfile(d, ['r' f '.mat']), 'file'))
                movefile(fullfile(d, ['r' f '.mat']), [resliced_img(1:end-4) '.mat']);
            end
        end
        
        if i==1
            reslicedImgs = resliced_img;
        else
            reslicedImgs = strvcat(reslicedImgs, resliced_img);
        end
    end
end

function reslice_name = new_resliced_name(overwrite, d, f, e)
    if(overwrite)
        reslice_name = fullfile(d, [f e]);
    else
        reslice_name = fullfile(d, ['r' f e]);
    end
end