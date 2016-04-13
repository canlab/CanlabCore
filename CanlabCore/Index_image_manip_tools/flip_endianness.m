function flip_endianness(imgs)
% Flips fMRI images (swap R/L)
%
% :Usage:
% ::
%
%     [Vimgs] = flip_endianness(imgs)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) <2006>  <CANLab>
%     Written by Tor? Documented by Marianne
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **imgs:**
%        Give this a list of fMRI images (paths, char str)
%
% :Outputs:
%
%   **Vimgs:**
%        Your images, flipped
%
% :Examples:
% ::
%    load EXPT  %study specific spm design file
%    imgs = EXPT.SNPM.P{1};
%    [Vimgs] = flip_endianness(imgs)
%
% :References:
%   N/A
%
% :See also:
%   - spm_vol
%
% ..
%    Programmers' notes:
%    This just capitalizes on spm tools
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