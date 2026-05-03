function voxel_data_series = plot_current_orthviews_coord(dat)
% plot_current_orthviews_coord Retrieve and plot the image data series at the current SPM orthviews crosshairs.
%
% Looks up the voxel under the current spm_orthviews crosshair, finds its
% row in the image_vector data, and plots the values across images
% (e.g., across subjects, contrasts, or time points). If the voxel is
% outside the in-mask voxel list, prints a notice in the figure.
%
% :Usage:
% ::
%
%     voxel_data_series = plot_current_orthviews_coord(dat)
%
% :Inputs:
%
%   **dat:**
%        An image_vector / fmri_data / statistic_image object with a
%        valid .volInfo and .xyzlist.
%
% :Outputs:
%
%   **voxel_data_series:**
%        A 1 x n_images vector of values from dat.dat at the orthviews
%        coordinate, or [] if the coordinate is outside the in-mask
%        voxels or no valid volInfo is available.
%
% :Examples:
% ::
%
%     orthviews(dat);
%     y = plot_current_orthviews_coord(dat);
%
% :See also:
%   - orthviews
%   - spm_orthviews
%   - mm2voxel

voxel_data_series = [];

if isempty(dat.volInfo) || ~isfield(dat.volInfo, 'mat') || isempty(dat.volInfo.mat)
    disp('Warning: could not find valid volInfo field in data object');
    return
end

XYZ = mm2voxel(spm_orthviews('Pos'), dat.volInfo.mat);

if isempty(XYZ)
    disp('Warning: could not retrieve coordinate from SPM orthviews');
    return
end

[isinmask, loc] = ismember(XYZ, dat.volInfo.xyzlist(~dat.removed_voxels, :), 'rows');

create_figure('image series');
if isinmask
    
    voxel_data_series = dat.dat(loc, :);
    plot(voxel_data_series, 'k.-');
    xlabel('Image series');
    ylabel('Value');
    
    %hold on; plot(dat2.dat(loc, :) + 100, 'r.-');
    
else
    
    text(0, 0, 'Vox is outside mask', 'FontSize', 24);
    set(gca, 'XLim', [-.5 1], 'YLim', [-1 1]);
    axis off
end

end
