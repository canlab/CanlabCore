cutoff = 50;
aspectr = 1;
original = 0;

Ds = smooth3(D);

if original
hiso = patch(isosurface(Ds,cutoff),...
    'FaceColor',[1,.75,.65],...
    'EdgeColor','none');

hcap = patch(isocaps(D,cutoff),...
    'FaceColor','interp',...
    'EdgeColor','none');
else
[hiso,hcap] = imagePatch(X,Y,Z,D,[nan coords(1) nan nan nan 40]);
end

view(45,30) 
axis tight 
daspect([1,1,aspectr])

myLight = lightangle(45,30); 
set(gcf,'Renderer','zbuffer'); lighting phong
isonormals(Ds,hiso)
set(hcap,'AmbientStrength',.6)
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)

% this is the key to avoiding dark surface head.
[az,el] = view
lightangle(az,el)