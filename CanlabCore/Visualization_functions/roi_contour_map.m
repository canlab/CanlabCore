function info = roi_contour_map(dat, varargin)
% Draw a pattern map of one slice (either saggital, axial, or coronal) that
% shows the most voxels, or the slice that you specify (e.g., x = #). 
% You can also draw outlines for the significant voxels from a statistical test. 
%
% :Usage:
% ::
%
%    info = roi_contour_map(dat, varargin)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2014 Choong-Wan (Wani) Woo
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
%   **dat:**
%        dat can be fmri_data, statistic_image (to mark significant 
%        voxels), and region objects. If your dat is the "region" 
%        object, please add 'cluster' as an optional input. 
%        You can display two pattern maps for the purpose of
%        comparison by putting additional columns of data in cell array. 
%
%        An example of displaying two pattern maps:
%        ::
%
%            dat{1} = region('img1.nii'); 
%            dat{2} = region('img2.nii'); )
%
% :Optional Inputs:
%
%   **'cluster':**
%        When the data is a region object, you need this option.
%
%   **'sig':**
%        This option outlines significant voxels. To use this
%        option, data in a format of statistic_image with a "sig" 
%        field should be given. 
%
%   **'colorbar':**
%        display colorbar under the plot. 
%
%   **'use_same_range':**
%        When you display two pattern maps, this option uses the
%        same color range for the two maps. 
%
%   **'surf':**
%        surface plot rather than voxel-by-voxel mapping.
%
%   **'xyz':**
%        When you want a specific view and slice, you can use this
%        option with 'coord'. (1:x - saggital view, 2:y - coronal 
%        view, 3:z - axial view)
%
%   **'coord':**
%        With 'xyz' option, this specifies the slice displayed. 
%
%   **'notfill':**
%        Default is to fill in the blank voxels using a black
%        color. With this option, you can color the blank voxels
%        with the white color. 
%
%   **'whole':**
%        Default is dividing the data into contiguous regions and
%        show only one region that has the most voxels. This option
%        akes this function not to divide into contiguous regions. 
%
%   **'colors' or 'color':**
%        you can specify your own colormap. 
%
%   **'contour':**
%        Not fully implemented yet. 
%
% :Outputs:
%
%   **info:**
%        information about the display with the following fields.
%
%   info.dat:
%     - [2x30 double] (xyz mesh)
%     - Z: [1x30 double] (z values)
%     - xyz: 3 (1:x-saggital, 2:y-coronal, 3:z-axial)
%     - xyz_coord: -2 (slice coordinate; in this case, z = 3)
%     - region_idx: 1 
%
%
% :Examples: you can also see the same example and output in 
% http://wagerlab.colorado.edu/wiki/doku.php/help/core/figure_gallery
% ::
%
%    mask{1} = 'dACC_hw_pattern_sl6mm.nii';
%    mask{2} = 'dACC_rf_pattern_sl6mm.nii';
%    for i = 1:2, cl{i} = region(mask{i}); end
%    info = roi_contour_map([cl{1} cl{2}], 'cluster', 'use_same_range', 'colorbar'); 
%
% ..
%    Copyright (C) 2014  Wani Woo
% ..

docolor = 0; % parsing optional inputs
dosig = 0;
docolorbar = 0; 
usesamerange = 0;
usecluster = 0;
dosurf = 0;
docontour = 0; % needs to be implemented.
doxyz = 0;
docoord = 0;
donotfill = 0;
dowhole = 0;
outline_color_pos = 'r';
outline_color_neg = 'b';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'colors', 'color'}
                docolor = 1; colors = varargin{i+1};
            case {'sig'}
                dosig = 1; 
            case {'sigcolor'}
                outline_color_pos = varargin{i+1};
                outline_color_neg = varargin{i+2};
            case {'colorbar'}
                docolorbar = 1; 
            case {'use_same_range'}
                usesamerange = 1; 
            case {'cluster'}
                usecluster = 1; 
            case {'surf'}
                dosurf = 1;
            case {'contour'}
                docontour = 1;
            case {'xyz'}
                doxyz = 1;
                xyz_idx = varargin{i+1};
            case {'coord'}
                docoord = 1;
                coord_idx = varargin{i+1};
            case {'notfill'}
                donotfill = 1;
            case {'whole'}
                dowhole = 1;
        end
    end
end

% data inputs
if ~usecluster
    if ~iscell(dat)
        dat_temp{1} = dat;
        dat = dat_temp;
    end
end

rnum  = numel(dat);

for jj = 1:rnum
    if ~usecluster
        if dosig
            if any(strcmp(fields(dat{jj}), 'sig'))
                datsig{jj} = dat{jj}; 
                datsig{jj}.dat = dat{jj}.sig; 
                datsig{jj}.dat = dat{jj}.dat.*~dat{jj}.sig*10^-6 + dat{jj}.sig;
                dat{jj}.sig = ones(size(dat{jj}.dat)); datsig{jj}.sig = ones(size(datsig{jj}.dat));
            else
                error('There is no sig data. Please check your data.');
            end
        else
            if any(strcmp(fields(dat{jj}), 'sig'))
                dat{jj}.sig = ones(size(dat{jj}.dat));
            end
        end
        
        ri = region(dat{jj});
        if dosig, risig = region(datsig{jj}); end
        if dowhole, ri = combine_region(ri); risig = combine_region(risig); end
    else
        ri = dat(jj);
    end
    clear max_vox idx idxidx;
    
    for i = 1:numel(ri)
        
        if ~docoord
            mm_c{i} = center_of_mass(ri(i).XYZmm, ones(1,size(ri(i).XYZmm, 2)));
        elseif docoord && ~doxyz
            warning('''coord'' option should be used with ''xyz''. This will use center of mass.');
            mm_c{i} = center_of_mass(ri(i).XYZmm, ones(1,size(ri(i).XYZmm, 2)));
        else
            mm_c{i} = center_of_mass(ri(i).XYZmm, ones(1,size(ri(i).XYZmm, 2)));
            [~, idx_coord_idx] = min(abs(ri(i).XYZmm(xyz_idx,:)-coord_idx));
            new_coord_idx = ri(i).XYZmm(xyz_idx,idx_coord_idx);
            mm_c{i}(xyz_idx) = new_coord_idx;
        end
        
        if ~doxyz
            [max_vox(i), idx(i)] = max([sum(ri(i).XYZmm(1,:) == mm_c{i}(1)), sum(ri(i).XYZmm(2,:) == mm_c{i}(2)), sum(ri(i).XYZmm(3,:) == mm_c{i}(3))]);
        else
            voxs = [sum(ri(i).XYZmm(1,:) == mm_c{i}(1)), sum(ri(i).XYZmm(2,:) == mm_c{i}(2)), sum(ri(i).XYZmm(3,:) == mm_c{i}(3))];
            idx(i) = xyz_idx;
            max_vox(i) = voxs(xyz_idx);
        end
        
    end
    
    [~, idxidx] = max(max_vox);
    
    xyz = 1:3;
    xyz(idx(idxidx)) = [];
    
    coord{jj}.dat = ri(idxidx).XYZmm(xyz,ri(idxidx).XYZmm(idx(idxidx),:) == mm_c{idxidx}(idx(idxidx)));
    coord{jj}.Z = ri(idxidx).Z(ri(idxidx).XYZmm(idx(idxidx),:) == mm_c{idxidx}(idx(idxidx)));
    if dosig
        coord{jj}.sig = risig(idxidx).Z(risig(idxidx).XYZmm(idx(idxidx),:) == mm_c{idxidx}(idx(idxidx)));
    end
    coord{jj}.xyz = idx(idxidx);
    coord{jj}.xyz_coord = mm_c{idxidx}(idx(idxidx));
    coord{jj}.region_idx = idxidx;
    
end

info = coord;

%% basic settingfor drawing

close all;
clear xyz;
figure;
if docolorbar
    set(gcf, 'color', 'w', 'position', [ 6   221   500*rnum   485]); % with colorbar
else
    set(gcf, 'color', 'w', 'position', [ 6   322   500*rnum   384]);
end

col = [0.0461    0.3833    0.5912
    0.2461    0.5833    0.7912
    0.4000    0.7608    0.6471
    0.6706    0.8667    0.6431
    0.9020    0.9608    0.5961
    1.0000    1.0000    0.7490
    0.9961    0.8784    0.5451
    0.9922    0.6824    0.3804
    0.9569    0.4275    0.2627
    0.8353    0.2431    0.3098
    0.6196    0.0039    0.2588];

colormap(col);
if docolor, colormap(colors); end

for jj = 1:rnum
    
    if dosig
        xyz = [coord{jj}.dat(1,:)', coord{jj}.dat(2,:)', coord{jj}.Z', sign(coord{jj}.Z') .* coord{jj}.sig'];
    else
        xyz = [coord{jj}.dat(1,:)', coord{jj}.dat(2,:)', coord{jj}.Z'];
    end
    
    [vX, vY] = meshgrid(unique(xyz(:,1)), unique(xyz(:,2)));
    
    vZ{jj} = NaN(size(vX));
    
    for i = 1:size(vX,2)
        for j = 1:size(vY,1)
            for k = 1:length(xyz)
                if isequal([vX(1,i) vY(j,1)], xyz(k,1:2))
                    vZ{jj}(j,i) = xyz(k,3);
                    if dosig
                        vsig{jj}(j,i) = xyz(k,4);
                    end
                end
            end
        end
    end
    
    vZ{jj} = flipud(vZ{jj});
    if dosig
        vsig{jj} = flipud(vsig{jj});
    end
    if docontour
        vZ{jj} = interp2(vZ{jj},5);
    end
    
    vox_size = vX(1,2)-vX(1,1);
    newX = vX/vox_size-vX(1,1)/vox_size+1;
    newY = vY/vox_size-vY(1,1)/vox_size+1;
    
    if docontour
        newX = repmat(1:size(vZ{jj},2), size(vZ{jj},1),1);
        newY = repmat((1:size(vZ{jj},1))', 1, size(vZ{jj},2));
    end
    
    %newX = (vX+2)/3-(vX(1,1)+2)/3+1;
    %newY = (vY+2)/3-(vY(1,1)+2)/3+1;
            
    xx{jj} = [newX(isnan(vZ{jj}))-.5 newX(isnan(vZ{jj}))-.5 newX(isnan(vZ{jj}))+.5 newX(isnan(vZ{jj}))+.5];
    yy{jj} = [newY(isnan(vZ{jj}))-.5 newY(isnan(vZ{jj}))+.5 newY(isnan(vZ{jj}))+.5 newY(isnan(vZ{jj}))-.5];
    
    if donotfill
        xx_outline{jj} = [newX(~isnan(vZ{jj}))-.5 newX(~isnan(vZ{jj}))-.5 newX(~isnan(vZ{jj}))+.5 newX(~isnan(vZ{jj}))+.5 newX(~isnan(vZ{jj}))-.5];
        yy_outline{jj} = [newY(~isnan(vZ{jj}))-.5 newY(~isnan(vZ{jj}))+.5 newY(~isnan(vZ{jj}))+.5 newY(~isnan(vZ{jj}))-.5 newY(~isnan(vZ{jj}))-.5];
    end
    
    if dosig
        x_out_pos{jj} = [newX(vsig{jj}==1)-.5 newX(vsig{jj}==1)-.5 newX(vsig{jj}==1)+.5 newX(vsig{jj}==1)+.5 newX(vsig{jj}==1)-.5];
        y_out_pos{jj} = [newY(vsig{jj}==1)-.5 newY(vsig{jj}==1)+.5 newY(vsig{jj}==1)+.5 newY(vsig{jj}==1)-.5 newY(vsig{jj}==1)-.5];
        x_out_neg{jj} = [newX(vsig{jj}==-1)-.5 newX(vsig{jj}==-1)-.5 newX(vsig{jj}==-1)+.5 newX(vsig{jj}==-1)+.5 newX(vsig{jj}==-1)-.5];
        y_out_neg{jj} = [newY(vsig{jj}==-1)-.5 newY(vsig{jj}==-1)+.5 newY(vsig{jj}==-1)+.5 newY(vsig{jj}==-1)-.5 newY(vsig{jj}==-1)-.5];
    end
    
    lim_max(jj) = max(max(vZ{jj}));
    lim_min(jj) = min(min(vZ{jj}));    
end

if usesamerange
    lim_max = max(lim_max); lim_min = min(lim_min);
end

for jj = 1:rnum
    
    subplot(1,rnum,jj);
    
    if dosurf % not work that well yet, but kind of work
        % surf(newX, newY, vZ{jj}, 'EdgeColor','k');
        surf(interp2(vZ{jj}, 5), 'EdgeColor','none');
        set(gca, 'view', [-21.5 70])
        % axis off;
        % camlight left; lighting phong;
        set(gca, 'linewidth', 1, 'xtick', [], 'ytick', [], 'ztick', [], 'fontsize', 15);
        if usesamerange
            set(gca, 'zlim', [lim_min lim_max]);
        end
    else
        if ~usesamerange
            if docontour
                imagesc(vZ{jj}, [lim_min(jj) lim_max(jj)]);
            else
                imagesc(vZ{jj}, [lim_min(jj) lim_max(jj)]);
                if ~donotfill
                    set(gca, 'linewidth', 5, 'xtick', [], 'ytick', []);
                else
                    set(gca, 'linewidth', 0.0001, 'xtick', [], 'ytick', [], 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w');
                end
            end
        else
            if docontour
                imagesc(vZ{jj}, [lim_min lim_max]);
                set(gca, 'linewidth', 0.0001, 'xtick', [], 'ytick', [], 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w');                    
            else
                imagesc(vZ{jj}, [lim_min lim_max]);
                if ~donotfill
                    set(gca, 'linewidth', 5, 'xtick', [], 'ytick', []);
                else
                    set(gca, 'linewidth', 0.0001, 'xtick', [], 'ytick', [], 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w');                    
                end
            end
        end
        
        if ~docontour
            y = size(vZ{jj},1); x = size(vZ{jj},2);
            if ~donotfill
                %for i = 1:x, line([i+.5,i+.5], [0 y+.5], 'color', [.2 .2 .2], 'linewidth', 1.5, 'linestyle', '-'); end
                %for i = 1:y, line([0 x+.5], [i+.5,i+.5], 'color', [.2 .2 .2], 'linewidth', 1.5, 'linestyle', '-'); end
                
                % black lines
                for i = 1:x, line([i+.5,i+.5], [0 y+.5], 'color', [.2 .2 .2], 'linewidth', 1, 'linestyle', '-'); end
                for i = 1:y, line([0 x+.5], [i+.5,i+.5], 'color', [.2 .2 .2], 'linewidth', 1, 'linestyle', '-'); end

                % white lines
%                 for i = 1:x, line([i+.5,i+.5], [0 y+.5], 'color', 'w', 'linewidth', 1, 'linestyle', '-'); end
%                 for i = 1:y, line([0 x+.5], [i+.5,i+.5], 'color', 'w', 'linewidth', 1, 'linestyle', '-'); end
            else
                for i = 1:x, line([i+.5,i+.5], [0 y+.5], 'color', repmat(.2, 1, 3), 'linewidth', 1.5, 'linestyle', '-'); end
                for i = 1:y, line([0 x+.5], [i+.5,i+.5], 'color', repmat(.2, 1, 3), 'linewidth', 1.5, 'linestyle', '-'); end
            end
        end
            
        hold on;
        if docontour
            for i = 1:size(xx{jj},1), fill(xx{jj}(i,:), yy{jj}(i,:), 'w'); end % this is wrong.
        else
            if ~donotfill
                % fill black
                for i = 1:size(xx{jj},1), fill(xx{jj}(i,:), yy{jj}(i,:), 'k'); end
                
                % fill white
%                 for i = 1:size(xx{jj},1), fill(xx{jj}(i,:), yy{jj}(i,:), 'w', 'edgecolor', 'w'); end
            else
                for i = 1:size(xx{jj},1), fill(xx{jj}(i,:), yy{jj}(i,:), 'w', 'edgecolor', 'w'); end
                xy_outline_idx{jj} = outline(xx_outline{jj}, yy_outline{jj});
                draw_outline(xx_outline{jj},yy_outline{jj},xy_outline_idx{jj},'k',5);
            end
        end
        
        if dosig
            xy_idx_pos{jj} = outline(x_out_pos{jj}, y_out_pos{jj});
            draw_outline(x_out_pos{jj},y_out_pos{jj},xy_idx_pos{jj},outline_color_pos,7)
            
            xy_idx_neg{jj} = outline(x_out_neg{jj}, y_out_neg{jj});
            draw_outline(x_out_neg{jj},y_out_neg{jj},xy_idx_neg{jj},outline_color_neg,7)
%             for i = 1:size(x_out{jj},1)
%                 for i2 = 1:4
%                     if xy_idx{jj}(i,i2) 
%                         line([x_out{jj}(i,i2) x_out{jj}(i,i2+1)], [y_out{jj}(i,i2) y_out{jj}(i,i2+1)], 'color', 'y', 'linewidth', 5); 
%                     end
%                 end
%             end
        end
    end
    if docolorbar
        hh = colorbar('southoutside');
        set(hh, 'fontSize', 25, 'lineWidth', 3);
    end
    
    info{jj}.vZ = vZ{jj};
end

end

% subfunctions
function xy_idx = outline(x, y)

xy_idx = true(size(x,1),4);

for i = 1:size(x,1)
    for ii = 1:4
        for j = i+1:size(x,1)
            for jj = 1:4
                if sum(abs([sort(x(i,(ii:ii+1))) sort(y(i,(ii:ii+1)))] - [sort(x(j,(jj:jj+1))), sort(y(j,(jj:jj+1)))]) < .0001) == 4
                    xy_idx(i,ii) = false; xy_idx(j,jj) = false;
                end
            end
        end
    end
end

end

function draw_outline(x,y,idx,color, linewidth)
for i = 1:size(x,1)
    for j = 1:4
        if idx(i,j)
            line([x(i,j) x(i,j+1)], [y(i,j) y(i,j+1)], 'color', color, 'linewidth', linewidth);
        end
    end
end

end

function r = combine_region(r)
rr = r(1);

if numel(r) > 1
    for i = 2:numel(r)
        rr.XYZ = [rr.XYZ r(i).XYZ];
        rr.XYZmm = [rr.XYZmm r(i).XYZmm];
        rr.val = [rr.val; r(i).val];
        rr.Z = [rr.Z r(i).Z];
    end
end

r = rr;

end
