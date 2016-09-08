function [p, mesh_struct] = isosurface(obj, mythresh)

[~, ~, mesh_struct] = reconstruct_image(obj);                % get volume data for slices
voldata = mesh_struct.voldata;

mycolor = [.5 .5 .5];
mysmoothbox = 3;
mygaussstd = 1;
%mythresh = 80;

V = smooth3(voldata, 'gaussian', mysmoothbox, mygaussstd);

Fvoldata = isosurface(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, V, mythresh);

p = patch('Faces',Fvoldata.faces,'Vertices',Fvoldata.vertices,'FaceColor',[.5 .5 .5], ...
    'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',.3,'SpecularExponent',200);

set(p, 'FaceColor', mycolor);

lighting gouraud
camlight right
camlight left
material dull

end