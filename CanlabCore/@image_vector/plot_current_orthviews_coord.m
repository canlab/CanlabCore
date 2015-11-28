function voxel_data_series = plot_current_orthviews_coord(dat)
% Retrieves and plots the image data series at the current crosshairs in spm_orthviews
%

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
