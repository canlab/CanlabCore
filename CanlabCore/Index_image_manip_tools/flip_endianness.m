function flip_endianness(imgs)
    if(iscellstr(imgs))
        imgs = char(imgs);
    end
    fprintf('Loading meta-info...\n');
    Vimgs = spm_vol(imgs);

    fprintf('Flipping... ');
    for i=1:length(Vimgs)
        status_string = sprintf('%d/%d', i, length(Vimgs));
        fprintf(status_string);

        %        if(isfield(Vimgs(i), 'dt'))
        %        else
        %            if(Vimgs(i).dim(4) >= 256)
        %                Vimgs(i).dim(4) = Vimgs(i).dim(4) / 256;
        %            else
        %                Vimgs(i).dim(4) = Vimgs(i).dim(4) * 256;
        %            end
        %            Vimgs(i).private.hdr.dime.datatype = Vimgs(i).dim(4);
        %        end
        %
        %        data = spm_read_vols(Vimgs(i));
        %        spm_write_vol(Vimgs(i), data);

        [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = spm_hread(Vimgs(i).fname);
        if(TYPE >= 256)
            TYPE = TYPE / 256;
        else
            TYPE = TYPE * 256;
        end
        spm_hwrite(Vimgs(i).fname,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);

        erase_string(status_string);
    end

    fprintf('\nFinished converting %d images.\n', length(Vimgs));
end