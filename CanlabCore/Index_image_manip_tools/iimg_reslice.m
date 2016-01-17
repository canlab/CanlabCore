function out = iimg_reslice(matchto, reslicethis, varargin)
% :Usage:
% ::
%
%     out = iimg_reslice(matchto, reslicethis, varargin)
%
%     out = iimg_reslice(matchto, reslicethis, 'write', 'outname', 'myimg.img')
%
% ..
%    tor wager
%    nov. 06
% ..

    % default flags
    % interp = 2 is trilinear
    flags = struct('interp', 2, 'vox', NaN, 'bb', NaN, 'wrap', [0 0 0], 'preserve', 0);
    % spm options for interpolation and wrapping
    d  = [flags.interp*[1 1 1]' flags.wrap(:)];

    write_vol = 0;
    if any(strcmp(varargin, 'write'))
        write_vol = 1;
    end

    % get output image name
    outname = 'resliced_image.img';
    argno = find(strcmp(varargin, 'outname'));
    if argno
        outname = varargin{argno + 1};
    end

    % get info of image to match to
    volInfo = iimg_read_img(matchto);


    % grid in space of target image
    [x1, x2] = ndgrid(1:volInfo.dim(1), 1:volInfo.dim(2));

    out = zeros(volInfo.dim(1:3));

    for i = 1	% for each image
        % get info for image to reslice
        volInfo2 = spm_vol(reslicethis);
        %[volInfo2, reslicethis] = iimg_read_img(reslicethis);

        if write_vol
            VO = volInfo2;
            VO.fname   = outname;
            
            switch(spm('Ver'))
                case 'SPM2'
                    VO.dim     = [volInfo.dim(1:3) volInfo2.dim(4)];
                case {'SPM5', 'SPM8'}
                    VO.dt = volInfo.dt;
                    VO.private.dat.fname = outname;
                    
                otherwise
                    error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
            end

            VO.mat     = volInfo.mat;
            VO.descrip = 'iimg_reslice - resliced - trilinear';
        end


        % final affine matrix
        % mat is affine mtx of image
        % Target\Object maps target to object
        affinemat = inv(volInfo.mat\volInfo2.mat);


        C = spm_bsplinc(volInfo2, d);

        % for each slice

        for x3 = 1:volInfo.dim(3)
            [tmp, y1, y2, y3] = getmask(affinemat, x1, x2, x3, volInfo2.dim(1:3), flags.wrap);
            out(:, :, x3)              = spm_bsplins(C, y1, y2, y3, d);
        end

        if write_vol
            spm_write_vol(VO, out);
        end
    end
end




function [Mask, y1, y2, y3] = getmask(M, x1, x2, x3, dim, wrp)
    tiny = 5e-2; % From spm_vol_utils.c
    y1   = M(1, 1)*x1+M(1, 2)*x2+(M(1, 3)*x3+M(1, 4));
    y2   = M(2, 1)*x1+M(2, 2)*x2+(M(2, 3)*x3+M(2, 4));
    y3   = M(3, 1)*x1+M(3, 2)*x2+(M(3, 3)*x3+M(3, 4));
    Mask = logical(ones(size(y1)));
    if ~wrp(1), Mask = Mask & (y1 >= (1-tiny) & y1 <= (dim(1)+tiny)); end;
    if ~wrp(2), Mask = Mask & (y2 >= (1-tiny) & y2 <= (dim(2)+tiny)); end;
    if ~wrp(3), Mask = Mask & (y3 >= (1-tiny) & y3 <= (dim(3)+tiny)); end;
end
