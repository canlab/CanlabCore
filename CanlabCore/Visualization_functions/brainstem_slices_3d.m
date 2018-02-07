function [p, caphandle, mesh_struct] = brainstem_slices_3d()

%mask = which('spm8_brainstem_pag.img');
mask = which('canlab_brainstem.img');
overlay = which('SPM8_colin27T1_seg.img');

dat = fmri_data(overlay);
dat = apply_mask(dat, mask);                     % mask with brainstem
dat = threshold(dat, [80 Inf], 'raw-between');  % threshold to clean up mask
orthviews(dat);

[voldata, ~, mesh_struct] = reconstruct_image(dat);  % get volume data for slices
voldata = mesh_struct.voldata;                  % rotated for compatibility with addbrain.m

hh = findobj(gcf, 'Type', 'text'); delete(hh);
spm_orthviews('XHairs', 'off');
% scn_export_papersetup(600);
% saveas(gcf, 'brainstem_orthviews.png');

%%

create_figure('3dslices');

% you can show all slices like this:
skip = 10;
thickness = 0;

% these determine which slices to show
z_mm = [-40 -11];

% get z = top slice
for i = 1:length(z_mm)
    xyz = mm2voxel([0 0 z_mm(i)], dat.volInfo.mat);
    my_slices(i) = xyz(3);
end

nslices = size(voldata, 3);

p = [];
caphandle = [];

for slices = my_slices %1:skip:nslices-thickness
    
    wh_slice = [slices slices + thickness];
    
    v = voldata(:, :, wh_slice);

    FVC = isocaps(mesh_struct.X(:,:,wh_slice), mesh_struct.Y(:,:,wh_slice), mesh_struct.Z(:,:,wh_slice), v, 80);
    
    if thickness > 0
        Fvoldata = isosurface(mesh_struct.X(:,:,wh_slice), mesh_struct.Y(:,:,wh_slice), mesh_struct.Z(:,:,wh_slice), v, mythresh);
        
        p(end+1) = patch('Faces',Fvoldata.faces,'Vertices',Fvoldata.vertices, ...
            'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);
    end
    
    caphandle(end+1) = patch(FVC, 'FaceColor', 'interp', 'EdgeColor', 'none');
    drawnow
    
end

axis image
view(129, 14)



%% brainstem surface

[p_tmp, mesh_struct] = isosurface(dat,'sd', 2, 'thresh', .1);
p = [p p_tmp];

set(gca, 'ZLim', [-50 0])

% p = addbrain('hires left');
% set(p, 'FaceAlpha', .2, 'FaceColor', [.5 .5 .5]);

end % function

